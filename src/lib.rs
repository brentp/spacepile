use anyhow::{anyhow, Error, Result};
use rust_htslib::bam::record::{Aux, Cigar, CigarString, CigarStringView};
use rust_htslib::bam::{Read, Reader, Record};
//use std::error::Error;
use std::ops::Index;
use std::string::String;

pub struct SpacePile {
    bam: Reader,
    grouping_tag: String,
    last_group: String,
    max_length: u16,
    record: Record,
    error: Option<Error>,
}

impl Iterator for SpacePile {
    type Item = Vec<Vec<u16>>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut parts = Vec::new();
        while let Some(r) = self.bam.read(&mut self.record) {
            let tag = self.record.aux(self.grouping_tag.as_bytes());
            let group = if let Ok(group) = tag {
                match group {
                    Aux::String(group) => group,
                    _ => {
                        self.error =
                            Some(anyhow!("expecting sring tag for {:?}", self.grouping_tag));
                        return None;
                    }
                }
            } else {
                return None;
            };
            if parts.len() > 0 && self.last_group != group {
                space(&parts, self.max_length);
                parts.clear();
                self.last_group = String::from(group);
            }
            let cig = self.record.cigar();
            parts.push(cig);
        }
        None
    }
}

pub fn iterate(bam_path: String, grouping_tag: String) -> Result<usize> {
    let mut bam = Reader::from_path(&bam_path)?;

    let mut sp = SpacePile {
        bam,
        grouping_tag: String::from(grouping_tag),
        last_group: String::from(";;;;;;;"),
        max_length: 152,
        record: Record::new(),
        error: None,
    };

    Ok(0)
}

enum Encoding {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
    End = -1,
    Space = -2,
}

pub fn space(cigars: &Vec<CigarStringView>, max_length: u16) -> Result<Vec<Vec<u16>>> {
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
    let max_len = std::cmp::min(
        max_length,
        offsets.iter().map(|c| c.len).max().unwrap() as u16,
    );

    let mut result: Vec<Vec<u16>> = offsets
        .iter()
        .map(|_| Vec::with_capacity(max_len as usize))
        .collect();
    let min_start = offsets.iter().map(|c| c.cigar.pos()).min().unwrap();
    // if some reads start later, prepend '$' to pad.
    offsets
        .iter()
        .enumerate()
        .filter(|(i, o)| o.cigar.pos() > min_start)
        .for_each(|(i, o)| {
            for _ in 0..(o.cigar.pos() - min_start) {
                result[i].push(Encoding::End as u16)
            }
        });

    let mut L = 0;
    while offsets.iter().any(|o| o.idx < o.len as u16) {
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
                if result[i].len() > L {
                    // continue
                } else if o.idx >= o.len {
                    result[i].push(Encoding::End as u16);
                } else if let Cigar::Ins(ilen) = o.cigar[o.op_i as usize] {
                    o.op_off += 1;
                    if o.op_off == ilen as u16 {
                        o.op_i += 1;
                        o.op_off = 0;
                    }
                    result[i].push(o.idx - o.dels);
                    o.idx += 1;
                } else {
                    result[i].push(Encoding::Space as u16)
                }
            });
        } else {
            // no insertion so we add the base or del at each position.
            offsets.iter_mut().enumerate().for_each(|(i, o)| {
                if result[i].len() > L {
                    // continue
                } else if o.idx >= o.len {
                    result[i].push(Encoding::End as u16);
                } else {
                    if let Cigar::Del(_) = o.cigar[o.op_i as usize] {
                        o.dels += 1;
                        result[i].push(Encoding::Space as u16);
                    } else {
                        result[i].push(o.idx - o.dels);
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

        L += 1;
    }

    Ok(result)
}

impl Index<usize> for CigTracker<'_> {
    type Output = u32;
    fn index(&self, i: usize) -> &u32 {
        &0
    }
}

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
        let obs = space(&cigs, 7).expect("oh no");
        let exp: Vec<Vec<u16>> = vec![
            vec![0, 1, 2, 3, 65535, 65535],
            vec![65535, 0, 1, 65534, 65534, 2],
        ];
        assert_eq!(obs, exp);

        let cigs = vec![
            make(vec![Cigar::Match(4)], 0), //                                 AC   TG
            make(vec![Cigar::Match(2), Cigar::Del(2), Cigar::Match(1)], 1), // $C   T C
            make(vec![Cigar::Match(2), Cigar::Ins(3), Cigar::Match(2)], 0), // ACGGGTG
        ];
        // AC   TG
        // $CG  TG
        // ACGGGTG
        let obs = space(&cigs, 7).expect("oh no");
        let exp = vec![
            vec![0, 1, 65534, 65534],
            vec![65535, 0, 65534, 65534],
            vec![0, 1, 2, 3, 4, 5],
        ];

        // TODO: handle insertion and deletion at same place. not sure this is correct.
        let exp = vec![
            vec![0, 1, 65534, 65534, 65534, 2, 3, 65535, 65535],
            vec![65535, 0, 65534, 65534, 65534, 1, 65534, 65534, 2],
            vec![0, 1, 2, 3, 4, 5, 6, 65535, 65535],
        ];
        assert_eq!(obs, exp);
    }
}
