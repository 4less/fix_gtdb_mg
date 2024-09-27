use std::cmp::Ordering;

use clap::Parser;
use fix_gtdb_mg::{common::Args, pairwise_leakage::{Leakage, NormGenes, TinyTaxID}};



fn main() {
    let args: Args = Args::parse();

    let leakage = Leakage::load(&args);
    let normalized_leakage = leakage.normalize_incoming();
    let mut vec = normalized_leakage.into_iter().collect::<Vec<(TinyTaxID, NormGenes)>>();
    vec.sort_by(|(a, ag), (b, bg)| ag.total().partial_cmp(&bg.total()).unwrap_or(Ordering::Equal));
    for (l, g) in vec {
        println!("{}\t{}", l, g)
    }
}