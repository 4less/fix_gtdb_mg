
use std::{collections::HashMap, fmt::{format, Display}, path::Path};

use clap::Parser;
use fix_gtdb_mg::common::{sam_file_iterator, Sam};


#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
#[command(arg_required_else_help(true))]
#[command(max_term_width = 120)] // term_width sets it fixed, max term_width can be smaller
pub struct Args {
    /// Input file (.sam|.sam.gz)
    #[arg(short = 'i', long = "input", default_value_t = String::default())]
    pub input: String,

    /// Mapq threshold (filter everything strictly below)
    #[arg(short = 'm', long = "min_mapq", default_value_t = 4)]
    pub min_mapq: u8,

    /// Minimum number of genes to keep.
    #[arg(short = 'g', long = "min_genes", default_value_t = 60)]
    pub min_genes: i32,

    /// Tolerate this many incoming leaked reads
    #[arg(short = 'g', long = "genes", default_value_t = 10)]
    pub max_leaked_reads: i32,
}


pub fn taxid_geneid(token: &str) -> Result<(usize, usize), Box<dyn std::error::Error>> {
    let mut parts = token.split('_');
    
    // Use next() to get the first two parts and check for their existence
    let first_part = parts.next().ok_or("Missing first part")?;
    let second_part = parts.next().ok_or("Missing second part")?;

    Ok((first_part.parse()?, second_part.parse()?))
}

type TaxID = usize;
type GeneID = usize;

#[derive(Default)]
pub struct Leaks {
    pub correct: usize,
    pub incoming: usize,
    pub outgoing: usize,
}

#[derive(Default)]
pub struct Species {
    pub id: TaxID,
    pub leaks: Vec<Option<Leaks>>,
}

impl Display for Species {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut s = format!("{}\t{}\t{}", self.id, self.num_good_genes(0), self.num_leaked_on_genes(0));

        self.leaks.iter().enumerate().for_each(|(idx, e)| {
            let tmp = match e {
                Some(e) => format!("\t{},c:{},i:{},o:{}", idx, e.correct, e.incoming, e.outgoing),
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
        if geneid >= self.leaks.len() {
            self.leaks.resize_with(geneid + 1, || None);
            self.leaks[geneid] = Some(Leaks::default());
        }
        let len = self.leaks.len();
        self.leaks[geneid].as_mut().expect(&format!("geneid {}, length leaks {}", geneid, len))
    }

    pub fn add_correct(&mut self, geneid: GeneID) {
        let leaks = self.get(geneid);
        leaks.correct += 1;
    }

    pub fn add_incorrect(&mut self, geneid: GeneID, incoming: bool) {
        let leaks = self.get(geneid);
        if incoming { leaks.incoming += 1 } else { leaks.outgoing += 1};
    }

    pub fn num_genes(&self) -> usize {
        self.leaks.iter().filter(|x| x.is_some()).count()
    }
    
    pub fn num_leaked_on_genes(&self, threshold: usize) -> usize {
        self.leaks.iter().
            filter(|x| { match x {
                Some(x) => x.incoming > threshold,
                None => false,
            }}).count()
    }    

    pub fn total_incoming_leaks(&self, threshold: usize) -> usize {
        self.leaks.iter().
            filter(|x| { match x {
                Some(x) => x.incoming > threshold,
                None => false,
            }}).fold(0, |acc, x| acc + x.as_ref().unwrap().incoming)
    }

    pub fn num_good_genes(&self, threshold: usize) -> usize {
        self.num_genes() - self.num_leaked_on_genes(threshold)
    }

}

pub struct GeneLeaks {
    species: HashMap<TaxID, Species>,
}

impl Default for GeneLeaks {
    fn default() -> Self {
        Self { species: Default::default() }
    }
}

impl GeneLeaks {
    pub fn count_correct(&mut self, species: TaxID, gene: GeneID) {
        let entry = self.species.entry(species).or_insert(Species::new(species));
        entry.add_correct(gene);
    }

    pub fn count_incorrect(&mut self, species: TaxID, gene: GeneID, incoming: bool) {
        let entry = self.species.entry(species).or_insert(Species::new(species));
        entry.add_incorrect(gene, incoming);
    }

    pub fn top_incoming(&self) -> Vec<(&TaxID, &Species)> {
        let mut result = self.species.iter().collect::<Vec<(&TaxID, &Species)>>();

        result.sort_by_key(|(id, s)| { (-(s.num_leaked_on_genes(0) as isize), -(s.total_incoming_leaks(0) as isize)) } );

        result
    }
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
            true => result.count_correct(query_tid, query_gid),
            false => {
                result.count_incorrect(ref_tid, ref_gid, true);
                result.count_incorrect(query_tid, query_gid, false);
            },
        }
    }

    
    result
}

fn main() {
    let args: Args = Args::parse();
    
    let leaks = get_gene_leaks(&args);

    for (id, s) in leaks.top_incoming().iter().rev() {
        eprintln!("{}", s);
    }
}

// Functions


