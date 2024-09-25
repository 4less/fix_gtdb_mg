use std::{collections::HashMap, fs::{read_to_string, File}, io::{self, BufRead}, path::Path};

// The output is wrapped in a Result to allow matching on errors.
// Returns an Iterator to the Reader of the lines of the file.
pub fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

pub fn get_labels_map(file: impl AsRef<Path>) -> (Vec<String>, HashMap<String, usize>) {
    let mut id2lab = Vec::new();
    let mut lab2id = HashMap::default();

    if let Ok(lines) = read_lines(file) {
        for line in lines {
            let line = line.expect("Corrupt file");

            let tokens = line.split("\t").collect::<Vec<&str>>();
            let id: usize = tokens[1].parse().unwrap();
            let lineage = tokens[3];

            let species = lineage.split(";").last().unwrap();

            if id >= id2lab.len() { 
                id2lab.resize_with(id+1, || String::default());
            };
            id2lab[id].push_str(species);
            lab2id.insert(species.to_string(), id);
        }
    }

    (id2lab, lab2id)
}