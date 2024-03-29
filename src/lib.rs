use anyhow::Result;
use ndarray::prelude::*;
use ndarray::ArrayViewMut2;
use numpy::{PyArray1, PyArray2};
use pyo3::prelude::*;
use pyo3::{
    exceptions, pymodule,
    types::{PyList, PyModule},
    PyResult, Python,
};
use rust_htslib::bam::record::{Cigar, CigarString, CigarStringView};

#[allow(dead_code)]
enum Encoding {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
    End = -1,
    Space = -2,
}

/// Given a set of Cigar+Positions, create an Array of indexes that can be used to align
/// the same set of variants.
pub fn space(cigars: &Vec<CigarStringView>, max_length: u32) -> Result<Array2<u32>> {
    let mut result = Array2::<u32>::zeros((cigars.len(), max_length as usize).f());

    let mut v = result.view_mut();
    space_fill(cigars, max_length, &mut v)?;
    Ok(result)
}

fn consumes_either(op: &Cigar, skip_soft: bool) -> bool {
    match op {
        Cigar::Pad(_) | Cigar::HardClip(_) => false,
        Cigar::SoftClip(_) => !skip_soft,
        _ => true,
    }
}

/// for each column in y, check if all rows are space except the label.
/// if so, remove the column and shift all columns to the right left by 1.
/// idxs tracks the lookup of a read column to the label.
pub fn remove_label_only_insertions_rs(
    y: &mut ArrayViewMut2<u32>,
    idxs: &mut ArrayViewMut1<u32>,
) -> Result<()> {
    let (rows, cols) = (y.shape()[0], y.shape()[1]);
    // now we get the final row which is the label row.
    let mut idx = 0;

    // iterate over each column in y and check if all rows are space except the label.
    let label_only_insertion_cols = y
        .axis_iter(Axis(1))
        .enumerate()
        .filter(|(i, col)| {
            let f = col
                .slice(s![..rows - 1])
                .iter()
                .all(|&v| v == Encoding::Space as u32);

            if !f {
                idxs[idx] = *i as u32;
                idx += 1;
            }
            f
        })
        .map(|(i, _)| i)
        .collect::<Vec<_>>();
    // since we increment after setting in the loop.
    // we check less than here instead of < cols - 1.
    while idx < cols {
        idxs[idx] = Encoding::End as u32;
        idx += 1;
    }

    // for each label_only_insertion, shift all the columns to the right left by 1.
    for col in label_only_insertion_cols {
        for mut row in y.slice_mut(s![..rows - 1, ..]).axis_iter_mut(Axis(0)) {
            // TODO: how to do this without copy here?
            let right = row.slice(s![col + 1..]).into_owned();
            right.assign_to(&mut row.slice_mut(s![col..cols - 1]));
            row[cols - 1] = Encoding::End as u32;
        }
    }

    Ok(())
}

pub fn space_fill(
    cigars: &Vec<CigarStringView>,
    max_length: u32,
    result: &mut ArrayViewMut2<u32>,
) -> Result<()> {
    result.fill(Encoding::End as u32);
    let mut offsets: Vec<CigTracker> = cigars
        .iter()
        .map(|c| {
            let mut t = CigTracker {
                cigar: c,
                query_idx: match c.first().unwrap() {
                    Cigar::SoftClip(n) | Cigar::Pad(n) => *n,
                    _ => 0,
                },
                dels: 0,
                op_i: (!consumes_either(c.first().unwrap(), true)) as u32, // we can skip the first hard-clip
                op_off: 0,
                query_len: c
                    .iter()
                    .map(|o| {
                        if consumes_either(o, false) {
                            o.len()
                        } else {
                            0
                        }
                    })
                    .sum(),
            };
            t.query_len -= match t.cigar.last().unwrap() {
                Cigar::SoftClip(n) => *n,
                _ => 0,
            };
            t
        })
        .collect();
    let min_start = offsets.iter().map(|c| c.cigar.pos()).min().unwrap();

    let mut sweep_col: usize = 0;
    while sweep_col < max_length as usize && offsets.iter().any(|o| o.query_idx < o.query_len) {
        let any_insertion = offsets.iter().any(|o| {
            if (o.op_i as usize) < o.cigar.iter().len() {
                matches!(o.cigar[o.op_i as usize], Cigar::Ins(_))
            } else {
                false
            }
        });
        if any_insertion {
            // found an insertion so we have to add a space in any read that did not have an
            // insertion. and add the base for the ones that did.
            offsets.iter_mut().enumerate().for_each(|(i, o)| {
                // this read started after. we already filled with Encoding::End above so here we
                // simply skip
                if o.cigar.pos() > min_start + sweep_col as i64 {
                    //if result[(i, sweep_col)] == Encoding::End as u32 {
                    // continue
                } else if o.query_idx - o.dels >= o.query_len || o.op_i >= o.cigar.len() as u32 {
                    //result[(i, sweep_col)] = Encoding::End as u32
                } else if let Cigar::Ins(ilen) = o.cigar[o.op_i as usize] {
                    o.op_off += 1;
                    if o.op_off == ilen {
                        o.op_i += 1;
                        o.op_off = 0;
                    }
                    result[(i, sweep_col)] = o.query_idx - o.dels;
                    o.query_idx += 1;
                } else {
                    result[(i, sweep_col)] = Encoding::Space as u32;
                }
            });
        } else {
            // no insertion so we add the base or del at each position.
            offsets.iter_mut().enumerate().for_each(|(i, o)| {
                if o.cigar.pos() > min_start + sweep_col as i64 {
                    // continue
                } else if o.query_idx - o.dels >= o.query_len || o.op_i >= o.cigar.len() as u32 {
                    //result[(i, sweep_col)] = Encoding::End as u32;
                } else {
                    if let Cigar::Del(_) = o.cigar[o.op_i as usize] {
                        o.dels += 1;
                        result[(i, sweep_col)] = Encoding::Space as u32;
                    } else {
                        // TODO: check the cigar op here? Only for Match() or Mismatch?
                        result[(i, sweep_col)] = o.query_idx - o.dels;
                    }
                    o.op_off += 1;
                    if o.op_off == o.cigar[o.op_i as usize].len() {
                        o.op_i += 1;
                        o.op_off = 0;
                    }
                    o.query_idx += 1;
                }
            });
        }

        sweep_col += 1;
    }

    Ok(())
}

/// Tracking for converting a group of cigars to a multi-alignment.
#[derive(Debug)]
struct CigTracker<'a> {
    cigar: &'a CigarStringView,
    query_len: u32,
    op_i: u32,
    /// tracks which op we are in.
    op_off: u32,
    query_idx: u32,
    /// total offset within the read.
    dels: u32,
}

/// create a rust-htslib cigar string from a python cigar string.
fn extract_cigar(cl: &PyAny, i: i64) -> CigarStringView {
    let cs: Vec<Cigar> = cl
        .downcast::<PyList>()
        .expect("expecting list as 3rd argument to 'space'")
        .iter()
        .map(|t| {
            let (op, len) = t.extract().unwrap();
            match op {
                0 => Cigar::Match(len),
                1 => Cigar::Ins(len),
                2 => Cigar::Del(len),
                3 => Cigar::RefSkip(len),
                4 => Cigar::SoftClip(len),
                5 => Cigar::HardClip(len),
                6 => Cigar::Pad(len),
                7 => Cigar::Equal(len),
                8 => Cigar::Diff(len),
                _ => panic!("Unexpected cigar operation"),
            }
        })
        .collect();
    CigarString(cs).into_view(i)
}

#[pymodule]
fn spacepile(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    /// fill array with spaced from py_cigs and positions. array.shape[0] == len(py_cigs) == len(py_positions) (or +1 if label is given)
    /// if a label position and cigs are provided, then arr.shape[0] == len(py_cigs) + 1.
    #[pyfunction]
    #[pyo3(name = "space")]
    fn rust_space<'py>(
        _py: Python<'py>,

        arr: &PyArray2<u32>,
        py_cigs: &PyList,
        py_positions: &PyList,
    ) -> PyResult<()> {
        // now get arr as a writeable, mutable array
        let mut y = unsafe { arr.as_array_mut() };

        if py_cigs.len() != py_positions.len() {
            return Err(exceptions::PyTypeError::new_err(
                "space: expecting py_cigs and py_positions of equal length",
            ));
        }
        if py_cigs.len() != arr.dims()[0] {
            return Err(exceptions::PyTypeError::new_err(
                "space: expecting py_cigs to have length equal to rows in arr",
            ));
        }

        let cigs: Vec<CigarStringView> = py_cigs
            .iter()
            .enumerate()
            .map(|(i, cl)| {
                extract_cigar(
                    cl,
                    py_positions
                        .get_item(i)
                        .expect("expecting element for position {i}")
                        .extract::<i64>()
                        .unwrap(),
                )
            })
            .collect();
        match space_fill(&cigs, arr.dims()[1] as u32, &mut y) {
            Ok(_) => Ok(()),
            Err(e) => Err(exceptions::PyTypeError::new_err(format!(
                "space: error filling array: {:?}",
                e
            ))),
        }
    }
    m.add_function(wrap_pyfunction!(rust_space, m)?)?;

    #[pyfunction]
    #[pyo3(name = "remove_label_only_insertions")]
    fn remove_label_only_insertions<'py>(
        _py: Python<'py>,
        spaced_idxs: &PyArray2<u32>,
        label_idxs: &PyArray1<u32>,
    ) -> PyResult<()> {
        // we remove any gaps in the label row that are not in the other rows.
        // a gap is a value of Encoding::Space as u32
        // we do this by iterating over the label row, and if we find a gap, we check if it is in any other row.
        // if it is not, we remove it. Now write the code for this.
        let mut y = unsafe { spaced_idxs.as_array_mut() };
        let mut label_is = unsafe { label_idxs.as_array_mut() };
        match remove_label_only_insertions_rs(&mut y, &mut label_is) {
            Ok(_) => Ok(()),
            Err(e) => Err(exceptions::PyTypeError::new_err(format!(
                "remove_label_only_insertions: error: {:?}",
                e
            ))),
        }
    }
    m.add_function(wrap_pyfunction!(remove_label_only_insertions, m)?)?;

    /// use `space` to translate  `sequences` (of sequence or BQs, for example) into `array`. `space` and `array` must be
    /// the same shape and `len(sequences) == space.shape[0]`
    #[pyfunction]
    #[pyo3(name = "translate")]
    fn translate<'py>(
        _py: Python<'py>,
        spaced_idxs: &PyArray2<u32>,
        sequences: &PyArray2<i16>,
        out: &PyArray2<i16>,
    ) -> PyResult<()> {
        if spaced_idxs.dims()[0] != sequences.dims()[0] {
            return Err(exceptions::PyTypeError::new_err(
                "translate: expecting spaced and sequences to have equal first dimension",
            ));
        }
        if spaced_idxs.dims()[0] != out.dims()[0] {
            return Err(exceptions::PyTypeError::new_err(
                "translate: expecting spaced and out to have equal first dimension",
            ));
        }
        let mut out_mut = unsafe { out.as_array_mut() };
        let spaced_raw = unsafe { spaced_idxs.as_array() };
        let sequences_raw = unsafe { sequences.as_array() };

        const UM1: u32 = u32::MAX - 1;

        for i in 0..spaced_idxs.dims()[0] {
            for j in 0..spaced_idxs.dims()[1] {
                let space_j = spaced_raw[(i, j)];
                out_mut[(i, j)] = match space_j {
                    u32::MAX => space_j as i16,
                    UM1 => space_j as i16,
                    _ => sequences_raw[(i, space_j as usize)],
                }
            }
        }
        Ok(())
    }

    m.add_function(wrap_pyfunction!(translate, m)?)
}

#[cfg(test)]
mod tests {
    use super::*;
    use pyo3::types::PyDict;
    use rust_htslib::bam::record::CigarString;

    fn make(cigs: Vec<Cigar>, pos: i64) -> CigarStringView {
        CigarStringView::new(CigarString(cigs), pos)
    }

    #[test]
    fn test_spacing() {
        let cigs = vec![
            make(vec![Cigar::Match(4)], 0), //                                 ACTG
            make(vec![Cigar::Match(2), Cigar::Del(2), Cigar::Match(1)], 1), // $CT  C
                                            //make(vec![Cigar::Match(2), Cigar::Ins(2), Cigar::Match(1)], 0),
        ];
        eprintln!("size:{:?}", std::mem::size_of::<CigTracker>());
        let obs = space(&cigs, 6).expect("oh no");
        let exp: Array2<u32> = array![
            [0, 1, 2, 3, u32::MAX, u32::MAX],
            [u32::MAX, 0, 1, u32::MAX - 1, u32::MAX - 1, 2],
        ];
        //                                                        [0, 1, 65534, 65534, 2, 65535]]
        assert_eq!(obs, exp);

        let cigs = vec![
            make(vec![Cigar::Match(4)], 0), //                                 AC   TG$$
            make(vec![Cigar::Match(2), Cigar::Del(2), Cigar::Match(1)], 1), // $C   T  C
            make(vec![Cigar::Match(2), Cigar::Ins(3), Cigar::Match(2)], 0), // ACGGGTG$$
        ];
        // AC   TG
        // $CG  TG
        // ACGGGTG
        let obs = space(&cigs, 9).expect("oh no");

        // TODO: handle insertion and deletion at same place. not sure this is correct.
        let exp: Array2<u32> = array![
            [
                0,
                1,
                u32::MAX - 1,
                u32::MAX - 1,
                u32::MAX - 1,
                2,
                3,
                u32::MAX,
                u32::MAX
            ],
            [
                u32::MAX,
                0,
                u32::MAX - 1,
                u32::MAX - 1,
                u32::MAX - 1,
                1,
                u32::MAX - 1,
                u32::MAX - 1,
                2
            ],
            [0, 1, 2, 3, 4, 5, 6, u32::MAX, u32::MAX],
        ];
        assert_eq!(obs, exp);
    }

    #[test]
    fn test_remove_label_only_insertions() {
        // create an array of shape (2, 5) with a gap in the middle
        let mut a = Array2::<u32>::zeros((3, 5));
        let mut idx = Array1::<u32>::zeros(5);
        // R1: AT-GC
        // R2: AA-GC
        // L : AAAGC
        a[(0, 2)] = Encoding::Space as u32;
        a[(1, 2)] = Encoding::Space as u32;
        remove_label_only_insertions_rs(&mut a.view_mut(), &mut idx.view_mut()).expect("oh no");
        let exp = array![
            [0, 0, 0, 0, Encoding::End as u32], // the 0's here are not important. real data
            [0, 0, 0, 0, Encoding::End as u32], // will have these as the base-encoded values.
            [0, 0, 0, 0, 0],
        ];
        assert_eq!(a, exp);

        // idx is a lookup from read_index to label index.
        assert_eq!(idx, array![0, 1, 3, 4, Encoding::End as u32]);
    }
    /*
        #[test]
        fn test_from_python() {
            pyo3::prepare_freethreaded_python();
            Python::with_gil(|py| {
                let locals = PyDict::new(py);
                let r = py.run(
                    r#"
    import numpy as np
    import spacepile
    import pysam

    a = np.zeros((2, 5), dtype=np.uint16)
    cigs = [[(int(pysam.CMATCH), 5)]]

    posns = [0]
    # NOTE we can send a slice
    spacepile.space(a[0:1], cigs, posns)
    assert np.array_equal(a[0].flatten(), np.arange(5, dtype=np.uint16))
            "#,
                    None,
                    Some(locals),
                );
                if r.is_err() {
                    eprintln!("error: {:?}", r.err());
                    assert!(false);
                } else {
                    assert!(r.is_ok());
                }
            });
        }

        #[test]
        fn test_translate_from_python() {
            pyo3::prepare_freethreaded_python();
            Python::with_gil(|py| {
                let locals = PyDict::new(py);
                let r = py.run(
                    r#"
    import numpy as np
    import spacepile
    import pysam

    a = np.zeros((2, 5), dtype=np.uint16)
    cigs = [[(pysam.CMATCH, 5)]]

    posns = [0]
    # NOTE we can send a slice
    spacepile.space(a[0:1], cigs, posns)
    mat = np.zeros((1, 5), dtype=np.int16)
    seqs = np.array([
        list('ACTTG')
        ], dtype='U1').view(np.int32).astype(np.int16)
    spacepile.translate(a[0:1], seqs, mat)
    assert "".join(chr(x) for x in mat[0]) == 'ACTTG', [chr(x) for x in mat[0]]
            "#,
                    None,
                    Some(locals),
                );
                if r.is_err() {
                    eprintln!("error: {:?}", r.err());
                    assert!(false);
                } else {
                    assert!(r.is_ok());
                }
            });
        }
        */
}
