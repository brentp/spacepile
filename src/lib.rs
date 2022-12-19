use anyhow::Result;
use ndarray::prelude::*;
use ndarray::Array2;
use rust_htslib::bam::record::{Cigar, CigarStringView};

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

    let mut result =
        Array2::<u16>::zeros((offsets.len(), max_length as usize).f()) + Encoding::End as u16;
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
                    result[(i, sweep_col)] = Encoding::End as u16
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
                    result[(i, sweep_col)] = Encoding::End as u16;
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

    Ok(result)
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

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::CigarString;

    fn make(cigs: Vec<Cigar>, pos: i64) -> CigarStringView {
        CigarStringView::new(CigarString(cigs), pos)
    }

    #[test]
    fn simple_spacing() {
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
}
