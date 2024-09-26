use std::{cmp::max, collections::HashMap, fmt::Display, path::Path};

use crate::common::{sam_file_iterator, sam_to_ids, Args, GeneID};



pub type TinyTaxID = u32;
pub type TinyGeneID = u32;

#[derive(Debug, Copy, Clone, PartialEq, PartialOrd, Eq, Hash)]
pub struct LeakagePair {
    pub from: TinyTaxID,
    pub to: TinyTaxID,
}

impl LeakagePair {
    pub fn from(origin: TinyTaxID, reference: TinyTaxID) -> Self {
        Self {
            from: origin,
            to: reference,
        }
    }
}

#[derive(Default)]
pub struct Genes {
    pub data: Vec::<isize>,
}

impl Display for Genes {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = itertools::join(self.data.iter().skip(1), "\t");
        write!(f, "{}\t{}", self.total(), s)
    }
}

impl Genes {
    const EMPTY: isize = -1;

    pub fn increment(&mut self, gene: GeneID) {
        if gene >= self.data.len() {
            self.data.resize_with(gene + 1, || Self::EMPTY);
        }

        self.data[gene] += 1
    }

    pub fn get(&self, gene: GeneID) -> Option<usize> {
        let count = self.data[gene];

        if count == Self::EMPTY {
            return None
        }
        return Some(count as usize)
    }

    pub fn total(&self) -> usize {
        self.data.iter().fold(0, |acc, x| acc + max(*x, 0)  as usize)
    }
}

#[derive(Default)]
pub struct Leakage {
    pub map: HashMap<LeakagePair, Genes>
}


impl Leakage {
    pub fn from_sam(args: &Args) -> Self {
        let mut iter = sam_file_iterator(&args.input).expect("Cannot open file");

        let mut res = Leakage::default();

        while let Some(Ok(sam)) = iter.next() {
            if !sam.is_aligned() || sam.mapq < args.min_mapq {continue};
            let fromto = sam_to_ids(&sam);
            
            let key = LeakagePair::from(fromto.query, fromto.reference);

            let mut entry = res.map.entry(key).or_default();

            if fromto.query_gene != fromto.reference_gene {
                eprintln!("Gene mismatch for Query Taxon: {} Gene: {} to Reference Taxon: {} Gene: {}", fromto.query, fromto.query_gene, fromto.reference, fromto.reference_gene);
            }

            entry.increment(fromto.reference_gene as GeneID);
        }
        res
    }
}