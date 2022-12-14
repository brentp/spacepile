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
        bam: bam,
        grouping_tag: String::from(grouping_tag),
        last_group: String::from(";;;;;;;"),
        max_length: 152,
        record: Record::new(),
        error: None,
    };

    Ok(0)
}

struct CigTracker<'a> {
    op_i: u16,
    /// tracks which op we are in.
    off: u16,
    /// total offset within the read.
    dels: u16,
    cigar: &'a CigarStringView,
    len: u32,
}

impl Index<usize> for CigTracker<'_> {
    type Output = u32;
    fn index(&self, i: usize) -> &u32 {
        // TODO:
        &0
    }
}

pub fn space(cigars: &Vec<CigarStringView>, max_length: u16) -> Result<Vec<Vec<u16>>> {
    let mut offsets: Vec<CigTracker> = cigars
        .iter()
        .map(|c| CigTracker {
            cigar: c,
            dels: 0,
            op_i: 0,
            off: 0,
            len: c.iter().map(|o| o.len()).sum(),
        })
        .collect();
    let max_len = [
        max_length,
        offsets.iter().map(|c| c.len).max().unwrap() as u16,
    ];

    while true {
        break;
    }

    Ok(vec![])
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make(cigs: Vec<Cigar>, pos: i64) -> CigarStringView {
        CigarStringView::new(CigarString(cigs), pos)
    }

    #[test]
    fn simple_spacing() {
        let a = make(vec![Cigar::Match(4)], 0);
        let b = make(vec![Cigar::Match(2), Cigar::Ins(3), Cigar::Match(2)], 0);
        let cigs = vec![a, b];
        space(&cigs, 7).expect("oh no");
    }
}
