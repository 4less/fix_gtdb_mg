use std::cmp::Ordering;

use clap::Parser;
use fix_gtdb_mg::{common::Args, pairwise_leakage::{Genes, Leakage, LeakagePair, NormGenes, TinyTaxID}};

fn main() {
    let args: Args = Args::parse();

    let leakage = Leakage::from_sam(&args);

    let mut vec = leakage.map.into_iter().collect::<Vec<(LeakagePair, Genes)>>();
    vec.sort_by_key(|l| (l.0.to, l.1.total()));
    for (l, g) in vec {
        println!("{}\t{}\t{}", l.from, l.to, g)
    }
}