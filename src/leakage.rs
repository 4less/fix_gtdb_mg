use std::{cmp::{max, min}, collections::HashMap, fmt::Display, fs::File, io::{BufRead, BufReader, Error}, iter::Map, path::Path};

use phylotree::tree::NodeId;

use crate::id_to_label::read_lines;


pub struct Leakage {
    pub from: NodeId,
    pub from_gene: NodeId,
    pub to: NodeId,
    pub to_gene: NodeId,
    pub correct: bool,
    pub mapq: usize,
}

impl Leakage {
    pub fn key(&self) -> (NodeId, NodeId) {
        return (min(self.from, self.to), max(self.from, self.to))
    }
}

#[derive(Default)]
pub struct LeakageCounter {
    pub total: usize,
    pub correct: usize,
    pub out_incorrect: usize,
    pub in_incorrect: usize,
}

impl Display for LeakageCounter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.total, 
            self.correct, self.correct as f64 / self.total as f64,
            self.out_incorrect, self.out_incorrect as f64 / self.total as f64,
            self.in_incorrect, self.in_incorrect as f64 / self.total as f64)
    }
}


// The output is wrapped in a Result to allow matching on errors.
// Returns an Iterator to the Reader of the lines of the file.
pub fn read_leakage_records<P>(filename: P) -> Map<std::io::Lines<BufReader<File>>, impl FnMut(Result<String, Error>) -> Leakage>
where P: AsRef<Path>, {
    let file = File::open(filename).unwrap();

    std::io::BufReader::new(file)
        .lines()
        .map(|line| {
            let line = line.expect("Corrupt file");

            let tokens = line.split("\t").collect::<Vec<&str>>();
            let from_tokens = tokens[1].split("_").collect::<Vec<&str>>();
            let to_tokens = tokens[2].split("_").collect::<Vec<&str>>();
            let correct = from_tokens == to_tokens;
            let mapq = tokens[4].parse().expect("Cannot parse bool");

            let from = from_tokens[0].parse().expect("Cannot parse bool");
            let from_gene = from_tokens[1].parse().expect("Cannot parse bool");
            let to = to_tokens[0].parse().expect("Cannot parse bool");
            let to_gene = to_tokens[1].parse().expect("Cannot parse bool");

            Leakage {
                from,
                from_gene,
                to,
                to_gene,
                correct,
                mapq,
            }
        }).into_iter()
}


pub fn read_leakage_file(path: impl AsRef<Path>) -> Vec<Leakage> {
    let mut result = Vec::default();
    let records = read_leakage_records(path);

    for record in records {
        result.push(record);
    }
    result
}

pub fn read_leakage_counter(path: impl AsRef<Path>) -> HashMap<NodeId, LeakageCounter> {
    let mut map = HashMap::new();

    let records = read_leakage_records(&path);
    for l in records {
        let from = map.entry(l.from).or_insert( LeakageCounter::default() );

        match l.correct {
            true => {
                from.correct += 1;
                from.total += 1;
            },
            false => {
                from.total += 1;
                from.out_incorrect += 1;
            }
        }
    }

    let records = read_leakage_records(path);
    for l in records {
        if l.correct {continue};
        let to = map.entry(l.to).or_insert( LeakageCounter::default() );
        to.in_incorrect += 1;
    }

    map
}

pub fn get_leakage_counter(leakage: &[Leakage]) -> HashMap<NodeId, LeakageCounter> {
    let mut map = HashMap::new();

    for l in leakage {
        let from = map.entry(l.from).or_insert( LeakageCounter::default() );

        match l.correct {
            true => {
                from.correct += 1;
                from.total += 1;
            },
            false => {
                from.total += 1;
                from.out_incorrect += 1;
            }
        }
    }

    for l in leakage {
        if l.correct {continue};
        let to: &mut LeakageCounter = map.entry(l.to).or_insert( LeakageCounter::default() );
        to.in_incorrect += 1;
    }

    map
}