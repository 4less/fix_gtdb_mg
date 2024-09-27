use std::{cmp::max, collections::HashMap, fmt::Display, path::Path};

use crate::{common::{sam_file_iterator, sam_to_ids, Args, GeneID}, utils::file_lines};



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

#[derive(Default)]
pub struct NormGenes {
    pub data: Vec::<f64>,
}

impl Display for Genes {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = itertools::join(self.data.iter().skip(1), "\t");
        write!(f, "{}\t{}", self.total(), s)
    }
}

impl Display for NormGenes {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = itertools::join(self.data.iter().skip(1), "\t");
        write!(f, "{}\t{}", self.total(), s)
    }
}

impl NormGenes {
    const EMPTY: f64 = -1.0;

    pub fn merge_normalized_from_counts(&mut self, other: &Genes, normalizer: &Genes) {
        for (gene, count) in other.data.iter().enumerate() {
            if *count == Genes::EMPTY { continue };



            if gene >= self.data.len() {
                self.data.resize_with(gene + 1, || Self::EMPTY);
            }
            if self.data[gene] == Self::EMPTY { self.data[gene] = 0.0 };

            let res = *count as f64 / normalizer.data[gene] as f64;

            assert!(res > 0.0);

            if res.is_nan() {
                eprintln!("Result: {}/{} = {}", *count as f64, normalizer.data[gene] as f64,  *count as f64 / normalizer.data[gene] as f64);
            }

            self.data[gene] += res;
        }
    }
    pub fn total(&self) -> f64 {
        let res = self.data.iter().fold(0.0, |acc, x| acc + if *x < 0.0 || *x == std::f64::NAN { 0.0 } else { *x }); //
        eprintln!("-- {} ... {} ... {:?}", res, res.is_nan(), self.data);
        assert!(res >= 0.0);

        res
    }
}

impl Genes {
    const EMPTY: isize = -1;

    pub fn increment(&mut self, gene: GeneID) {
        if gene >= self.data.len() {
            self.data.resize_with(gene + 1, || Self::EMPTY);
        }
        if self.data[gene] == Self::EMPTY { self.data[gene] = 0 };
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

    pub fn merge_from(&mut self, other: &Self) {
        for (gene, count) in other.data.iter().enumerate() {
            if *count == Self::EMPTY { continue };
            assert!(*count > 0);

            if gene >= self.data.len() {
                self.data.resize_with(gene + 1, || Self::EMPTY);
            }
            if self.data[gene] == Self::EMPTY { self.data[gene] = 0 };

            self.data[gene] += count;
            assert!(self.data[gene] > 0);
        }
    }

    pub fn from_slice(slice: &[isize]) -> Self {
        Self {
            data: Vec::from(slice),
        }
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

        while let Some(sam_res) = iter.next() {
            let sam = sam_res.expect("Invalid sam");
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
    
    pub fn load(args: &Args) -> Self {
        let mut result = Self::default();
        let mut tokens = Vec::<String>::new();
        let mut iter = file_lines(&args.input).expect("Unable to construct line iterator over file");
        while let Some(Ok(line)) = iter.next() {
            let mut tokens = line.split("\t").map(|x| x.parse().expect("Cannot parse string into i32")).collect::<Vec<isize>>();

            // eprintln!("{:?}", tokens);

            let from = tokens[0] as u32;
            let to = tokens[1] as u32;

            eprintln!("{:?}", &tokens[3..]);
            assert!(&tokens[3..].iter().all(|x| *x != 0));

            let key = LeakagePair::from(from, to);
            result.map.insert(key, Genes::from_slice(&tokens[3..]));
        }

        result
    }

    pub fn total_outgoing(&self) -> HashMap<TinyTaxID, Genes> {
        let mut result = HashMap::default();

        for (pair, genes) in &self.map {
            let from = pair.from;
            let entry: &mut Genes = result.entry(from).or_default();
            entry.merge_from(genes);
        }

        result
    }

    pub fn normalize_incoming(&self) -> HashMap<TinyTaxID, NormGenes>{
        let total_out = self.total_outgoing();
        let mut result = HashMap::default();

        for (pair, genes) in &self.map {
            let to: u32 = pair.to;
            let normalizer = &total_out[&pair.from];

            let entry: &mut NormGenes = result.entry(to).or_default();
            entry.merge_normalized_from_counts(genes, normalizer);

        }

        result
    }
}