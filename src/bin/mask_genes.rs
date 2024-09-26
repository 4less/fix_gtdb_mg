
use std::{collections::HashMap, fmt::Display, process::exit};

use clap::Parser;
use fix_gtdb_mg::common::{sam_file_iterator, taxid_geneid, Args, GeneID, TaxID};






#[derive(Default)]
pub struct Leaks {
    pub correct: f64,
    pub incoming: f64,
    pub outgoing: f64,
}

#[derive(Default)]
pub struct Species {
    pub id: TaxID,
    pub leaks: Vec<Option<Leaks>>,
}

impl Display for Species {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut s = String::default();

        s.push_str(&format!("{}\t{}\t{}\tcorrect", self.id, self.num_good_genes(0.0), self.num_leaked_on_genes(0.0)));
        self.leaks.iter().skip(1).for_each(|e: &Option<Leaks>| {
            let tmp = match e {
                Some(e) => { format!("\t{}", e.correct) },
                None => "\tNone".to_string(),
            };
            s.push_str(&tmp)
        });        
        s.push_str(&format!("\n{}\t{}\t{}\tincoming", self.id, self.num_good_genes(0.0), self.num_leaked_on_genes(0.0)));
        self.leaks.iter().skip(1).for_each(|e| {
            let tmp = match e {
                Some(e) => { format!("\t{}", e.incoming) },
                None => "\tNone".to_string(),
            };
            s.push_str(&tmp)
        });        
        s.push_str(&format!("\n{}\t{}\t{}\toutgoing", self.id, self.num_good_genes(0.0), self.num_leaked_on_genes(0.0)));
        self.leaks.iter().skip(1).for_each(|e| {
            let tmp = match e {
                Some(e) => { format!("\t{}", e.outgoing) },
                None => "\tNone".to_string(),
            };
            s.push_str(&tmp)
        });

        


        write!(f, "{}", s)
    }
}

impl Species {
    pub fn new(taxid: TaxID) -> Self {
        Self {
            id: taxid,
            leaks: Vec::new(),
        }
    }

    pub fn get(&mut self, geneid: GeneID) -> &mut Leaks {
        if geneid >= self.leaks.len() || self.leaks[geneid].is_none() {
            self.leaks.resize_with(geneid + 1, || None);
            self.leaks[geneid] = Some(Leaks::default());
        }
        let len = self.leaks.len();
        self.leaks[geneid].as_mut().expect(&format!("geneid {}, length leaks {}", geneid, len))
    }

    pub fn add_correct(&mut self, geneid: GeneID, increment: f64) {
        let leaks = self.get(geneid);
        leaks.correct += increment;
    }

    pub fn add_incorrect(&mut self, geneid: GeneID, incoming: bool, increment: f64) {
        let leaks = self.get(geneid);
        if incoming { leaks.incoming += increment } else { leaks.outgoing += increment};
    }

    pub fn num_genes(&self) -> usize {
        self.leaks.iter().filter(|x| x.is_some()).count()
    }
    
    pub fn num_leaked_on_genes(&self, threshold: f64) -> usize {
        self.leaks.iter().
            filter(|x| { match x {
                Some(x) => x.incoming > threshold,
                None => false,
            }}).count()
    }    

    pub fn total_incoming_leaks(&self, threshold: f64) -> f64 {
        self.leaks.iter().
            filter(|x| { match x {
                Some(x) => x.incoming > threshold,
                None => false,
            }}).fold(0.0, |acc, x| acc + x.as_ref().unwrap().incoming)
    }

    pub fn num_good_genes(&self, threshold: f64) -> usize {
        self.num_genes() - self.num_leaked_on_genes(threshold)
    }

}

pub struct GeneLeaks {
    species: HashMap<TaxID, Species>,
}

// type DirectionalLeakageKey = (TaxID, TaxID);

// pub struct DirectionalLeak {

// }

impl Default for GeneLeaks {
    fn default() -> Self {
        Self { species: Default::default() }
    }
}


impl GeneLeaks {
    pub fn count_correct(&mut self, species: TaxID, gene: GeneID, increment: f64) {
        let entry = self.species.entry(species).or_insert(Species::new(species));
        entry.add_correct(gene, increment);
    }

    pub fn count_incorrect(&mut self, species: TaxID, gene: GeneID, incoming: bool, increment: f64) {
        let entry = self.species.entry(species).or_insert(Species::new(species));
        entry.add_incorrect(gene, incoming, increment);
    }

    pub fn top_incoming(&self) -> Vec<(&TaxID, &Species)> {
        let mut result = self.species.iter().collect::<Vec<(&TaxID, &Species)>>();

        result.sort_by_key(|(_id, s)| { (-(s.num_leaked_on_genes(0.0) as isize), -(s.total_incoming_leaks(0.0) as isize)) } );

        result
    }
}



pub fn get_species_total(args: &Args) -> HashMap<TaxID, Vec<Option<usize>>> {
    let mut result = HashMap::default();

    
    let mut iter = sam_file_iterator(&args.input).expect("Cannot open file");


    while let Some(Ok(sam)) = iter.next() {
        // eprintln!("{:?}", sam);
        if !sam.is_aligned() || sam.mapq < args.min_mapq {continue};

        let (query_tid, query_gid) = taxid_geneid(&sam.qname).expect("Reference not parseable");

        let entry: &mut Vec<Option<usize>> = result.entry(query_tid).or_insert(Vec::default());
        if query_gid >= entry.len() || entry[query_gid].is_none() {
            entry.resize_with(query_gid + 1, || None);
            entry[query_gid] = Some(0);
        }
        *entry[query_gid].as_mut().unwrap() += 1;
    }

    result
}

pub fn get_normalized_gene_leaks(args: &Args, total_counts: &HashMap<TaxID, Vec<Option<usize>>>) -> GeneLeaks {
    let mut result = GeneLeaks::default();

    let mut iter = sam_file_iterator(&args.input).expect("Cannot open file");


    while let Some(Ok(sam)) = iter.next() {
        // eprintln!("{:?}", sam);
        if !sam.is_aligned() || sam.mapq < args.min_mapq {continue};

        let (query_tid, query_gid) = taxid_geneid(&sam.qname).expect("Reference not parseable");
        let (ref_tid, ref_gid) = taxid_geneid(&sam.rname).expect("Reference not parseable");
        let correct = query_tid == ref_tid && query_gid == ref_gid;

        let qt = total_counts.get(&query_tid);
        let query_total = match qt {
            Some(qt) => qt[query_gid].unwrap(),
            None => 0,
        };

        let rt = total_counts.get(&ref_tid);
        let ref_total = match qt {
            Some(rt) => rt[ref_gid].unwrap(),
            None => 0,
        };

        match correct {
            true => result.count_correct(query_tid, query_gid, 1.0 / query_total as f64),
            false => {
                result.count_incorrect(ref_tid, ref_gid, true, 1.0 / query_total as f64);
                result.count_incorrect(query_tid, query_gid, false, 1.0 / ref_total as f64);
            },
        }
    }

    
    result
}


pub fn get_gene_leaks(args: &Args) -> GeneLeaks {
    let mut result = GeneLeaks::default();

    let mut iter = sam_file_iterator(&args.input).expect("Cannot open file");


    while let Some(Ok(sam)) = iter.next() {
        // eprintln!("{:?}", sam);
        if !sam.is_aligned() || sam.mapq < args.min_mapq {continue};

        let (query_tid, query_gid) = taxid_geneid(&sam.qname).expect("Reference not parseable");
        let (ref_tid, ref_gid) = taxid_geneid(&sam.rname).expect("Reference not parseable");
        let correct = query_tid == ref_tid && query_gid == ref_gid;

        match correct {
            true => result.count_correct(query_tid, query_gid, 1.0),
            false => {
                result.count_incorrect(ref_tid, ref_gid, true, 1.0);
                result.count_incorrect(query_tid, query_gid, false, 1.0);
            },
        }
    }

    
    result
}


fn main() {
    let args: Args = Args::parse();
    
    let total = get_species_total(&args);

    let leaks = get_normalized_gene_leaks(&args, &total);

    eprintln!("{:?}", total);

    for (_id, s) in leaks.top_incoming().iter().rev() {
        println!("{}", s);
        
    }
}

// Functions


