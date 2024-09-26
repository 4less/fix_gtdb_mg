use clap::Parser;
use fix_gtdb_mg::{common::Args, pairwise_leakage::{Genes, Leakage, LeakagePair, NormGenes, TinyTaxID}};

fn main() {
    let args: Args = Args::parse();

    let leakage = Leakage::from_sam(&args);

    // let mut vec = leakage.map.into_iter().collect::<Vec<(LeakagePair, Genes)>>();
    // vec.sort_by_key(|l| (l.0.to, l.1.total()));
    // for (l, g) in vec {
    //     println!("{}\t{}\t{}", l.from, l.to, g)
    // }

    let normalized_leakage = leakage.normalize_incoming();
    let mut vec = normalized_leakage.into_iter().collect::<Vec<(TinyTaxID, NormGenes)>>();
    vec.sort_by(|(a, ag), (b, bg)| ag.total().partial_cmp(&bg.total()).unwrap());
    for (l, g) in vec {
        println!("{}\t{}", l, g)
    }
}