use std::{collections::HashMap, ffi::OsStr, fs::File, io::{BufRead, BufReader}, path::Path};

use clap::{command, Parser};
use flate2::bufread::GzDecoder;
use thiserror::Error;

use crate::pairwise_leakage::{TinyGeneID, TinyTaxID};

pub type TaxID = usize;
pub type GeneID = usize;


pub fn taxid_geneid(token: &str) -> Result<(usize, usize), Box<dyn std::error::Error>> {
    let mut parts = token.split('_');
    
    // Use next() to get the first two parts and check for their existence
    let first_part = parts.next().ok_or("Missing first part")?;
    let second_part = parts.next().ok_or("Missing second part")?;

    Ok((first_part.parse()?, second_part.parse()?))
}


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


/// Custom error type to handle different kinds of errors
#[derive(Debug, Error)]
pub enum SamFileError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("Gzip decode error: {0}")]
    Gzip(#[from] flate2::DecompressError),
}

/// Define a struct to represent a line in the SAM file
#[derive(Debug)]
pub struct Sam {
    pub qname: String,
    pub flag: u16,
    pub rname: String,
    pub pos: u32,
    pub mapq: u8,
    pub cigar: String,
    pub rnext: String,
    pub pnext: u32,
    pub tlen: i32,
    pub seq: String,
    pub qual: String,
}

impl Sam {
    /// Create a new Sam struct from a SAM file line
    pub fn from_line(line: &str) -> Result<Self, String> {
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 11 {
            return Err(format!("Invalid SAM line: {}", line));
        }

        Ok(Sam {
            qname: fields[0].to_string(),
            flag: fields[1].parse().map_err(|_| "Invalid flag")?,
            rname: fields[2].to_string(),
            pos: fields[3].parse().map_err(|_| "Invalid position")?,
            mapq: fields[4].parse().map_err(|_| "Invalid mapping quality")?,
            cigar: fields[5].to_string(),
            rnext: fields[6].to_string(),
            pnext: fields[7].parse().map_err(|_| "Invalid pnext")?,
            tlen: fields[8].parse().map_err(|_| "Invalid template length")?,
            seq: fields[9].to_string(),
            qual: fields[10].to_string(),
        })
    }

    pub fn is_aligned(&self) -> bool {
        return self.rname != "*";
    }
}

// type SamFileIterator = Result<impl Iterator<Item = Result<Sam, std::io::Error>>, SamFileError>;

/// A function that returns an iterator over Sam structs from a SAM file
pub fn sam_file_iterator<P: AsRef<Path>>(filename: P) ->  Result<impl Iterator<Item = Result<Sam, std::io::Error>>, SamFileError> {
    // Open the file
    let file = File::open(&filename)?;

    // Determine if the file is gzipped based on extension
    let reader: Box<dyn BufRead> = if filename.as_ref().extension() == Some(OsStr::new("gz")) {
        let gz_decoder = GzDecoder::new(BufReader::new(file));
        Box::new(BufReader::new(gz_decoder))
    } else {
        Box::new(BufReader::new(file))
    };

    // Create an iterator that processes each line into a Sam struct
    Ok(reader.lines().filter_map(|line_result| {
        match line_result {
            Ok(line) => {
                if line.starts_with('@') {
                    // Skip header lines starting with '@'
                    None
                } else {
                    // Parse the line into a Sam struct
                    match Sam::from_line(&line) {
                        Ok(sam) => Some(Ok(sam)),
                        Err(e) => Some(Err(std::io::Error::new(std::io::ErrorKind::InvalidData, e))),
                    }
                }
            }
            Err(e) => Some(Err(e)), // Propagate the I/O error
        }
    }))
}


/// Function to read a TSV file into a HashMap
pub fn read_tsv_to_hashmap<P: AsRef<Path>>(filename: P) -> std::io::Result<HashMap<String, String>> {
    // Open the file
    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    // Create an empty HashMap
    let mut hashmap = HashMap::new();

    // Read the file line by line
    for line_result in reader.lines() {
        let line = line_result?; // Handle any I/O error
        
        // Split the line by tab characters
        let mut columns = line.splitn(2, '\t');

        // Collect the first and second columns
        if let (Some(key), Some(value)) = (columns.next(), columns.next()) {
            hashmap.insert(key.to_string(), value.to_string());
        } else {
            eprintln!("Skipping invalid line: {}", line); // Handle cases where the line does not have exactly two columns
        }
    }

    Ok(hashmap)
}

pub struct FromTo {
    pub query: TinyTaxID,
    pub reference: TinyTaxID,
    pub query_gene: TinyGeneID,
    pub reference_gene: TinyGeneID,
}

pub fn sam_to_ids(sam: &Sam) -> FromTo {
    let (query_tid, query_gid) = taxid_geneid(&sam.qname).expect("Reference not parseable");
    let (ref_tid, ref_gid) = taxid_geneid(&sam.rname).expect("Reference not parseable");
    let correct = query_tid == ref_tid && query_gid == ref_gid;

    FromTo {
        query: query_tid as TinyTaxID,
        reference: ref_tid as TinyTaxID,
        query_gene: query_gid as TinyGeneID,
        reference_gene: ref_gid as TinyGeneID,
    }
}