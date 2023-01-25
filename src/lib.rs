use anyhow::Result;
use ndarray::prelude::*;
use ndarray::ArrayViewMut2;
use numpy::PyArray2;
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
pub fn space(cigars: &Vec<CigarStringView>, max_length: u16) -> Result<Array2<u16>> {
    let mut result =
        Array2::<u16>::zeros((cigars.len(), max_length as usize).f()) + Encoding::End as u16;

    let mut v = result.view_mut();
    space_fill(cigars, max_length, &mut v);
    Ok(result)
}

pub fn space_fill(
    cigars: &Vec<CigarStringView>,
    max_length: u16,
    result: &mut ArrayViewMut2<u16>,
) -> Result<()> {
    let mut offsets: Vec<CigTracker> = cigars
        .iter()
        .map(|c| CigTracker {
            cigar: c,
            idx: 0,
            dels: 0,
            op_i: 0,
            op_off: 0,
            len: c.iter().map(|o| o.len() as u16).sum(),
        })
        .collect();
    let min_start = offsets.iter().map(|c| c.cigar.pos()).min().unwrap();

    let mut sweep_col: usize = 0;
    while sweep_col < max_length as usize && offsets.iter().any(|o| o.idx < o.len as u16) {
        let any_insertion = offsets.iter().any(|o| {
            if (o.op_i as usize) < o.cigar.iter().len() {
                if let Cigar::Ins(_) = o.cigar[o.op_i as usize] {
                    true
                } else {
                    false
                }
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
                    //if result[(i, sweep_col)] == Encoding::End as u16 {
                    // continue
                } else if o.idx >= o.len {
                    //result[(i, sweep_col)] = Encoding::End as u16
                } else if let Cigar::Ins(ilen) = o.cigar[o.op_i as usize] {
                    o.op_off += 1;
                    if o.op_off == ilen as u16 {
                        o.op_i += 1;
                        o.op_off = 0;
                    }
                    result[(i, sweep_col)] = o.idx - o.dels;
                    o.idx += 1;
                } else {
                    result[(i, sweep_col)] = Encoding::Space as u16;
                }
            });
        } else {
            // no insertion so we add the base or del at each position.
            offsets.iter_mut().enumerate().for_each(|(i, o)| {
                if o.cigar.pos() > min_start + sweep_col as i64 {
                    // continue
                } else if o.idx >= o.len {
                    //result[(i, sweep_col)] = Encoding::End as u16;
                } else {
                    if let Cigar::Del(_) = o.cigar[o.op_i as usize] {
                        o.dels += 1;
                        result[(i, sweep_col)] = Encoding::Space as u16;
                    } else {
                        result[(i, sweep_col)] = o.idx - o.dels;
                    }
                    o.op_off += 1;
                    if o.op_off == o.cigar[o.op_i as usize].len() as u16 {
                        o.op_i += 1;
                        o.op_off = 0;
                    }
                    o.idx += 1;
                }
            });
        }

        sweep_col += 1;
    }

    Ok(())
}

/// Tracking for converting a group of cigars to a multi-alignment.
struct CigTracker<'a> {
    cigar: &'a CigarStringView,
    len: u16,
    op_i: u16,
    /// tracks which op we are in.
    op_off: u16,
    idx: u16,
    /// total offset within the read.
    dels: u16,
}

#[pymodule]
fn spacepile(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    /// space(array, py_cigs, py_positions, /)
    /// --
    ///
    /// fill array with spaced from py_cigs and positions. array.shape[0] == len(py_cigs) == len(py_positions)
    #[pyfn(m)]
    #[pyo3(name = "space")]
    fn rust_space<'py>(
        _py: Python<'py>,
        arr: &PyArray2<u16>,
        py_cigs: &PyList,
        py_positions: &PyList,
    ) -> PyResult<()> {
        //let mut y = x.readwrite();
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
                CigarString(cs).into_view(
                    py_positions
                        .get_item(i)
                        .expect("expecting element for position {i}")
                        .extract::<i64>()
                        .unwrap(),
                )
            })
            .collect();
        //let mut cigs = vec![];
        _ = space_fill(&cigs, 128, &mut y);
        Ok(())
    }

    /// translate(space, array, /)
    /// --
    ///
    /// use `space` to translate  `sequences` (of sequence or BQs, for example) into `array`. `space` and `array` must be
    /// the same shape and `len(sequences) == space.shape[0]`
    #[pyfn(m)]
    #[pyo3(name = "translate")]
    fn translate<'py>(
        _py: Python<'py>,
        spaced_idxs: &PyArray2<u16>,
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

        for i in 0..spaced_idxs.dims()[0] {
            for j in 0..spaced_idxs.dims()[1] {
                let space_j = spaced_raw[(i, j)];
                let si = sequences_raw[(i, space_j as usize)];
                out_mut[(i, j)] = si;
            }
        }
        Ok(())
    }

    Ok(())
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
        let exp: Array2<u16> = array![[0, 1, 2, 3, 65535, 65535], [65535, 0, 1, 65534, 65534, 2],];
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
        let exp: Array2<u16> = array![
            [0, 1, 65534, 65534, 65534, 2, 3, 65535, 65535],
            [65535, 0, 65534, 65534, 65534, 1, 65534, 65534, 2],
            [0, 1, 2, 3, 4, 5, 6, 65535, 65535],
        ];
        assert_eq!(obs, exp);
    }

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
}
