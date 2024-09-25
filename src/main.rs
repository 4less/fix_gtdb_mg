use std::{collections::{HashMap, HashSet}, env, fs::File, io::{BufReader, BufWriter, Write}, path::Path, process};

use id_to_label::get_labels_map;
use leakage::{get_leakage_counter, read_leakage_counter, read_leakage_file, Leakage};
use phylotree::tree::{Edge, Node, NodeId, Tree, TreeError};

pub mod id_to_label;
pub mod leakage;

fn parse_label(label: &str) -> String {
    String::default()
}

trait TreeHelper {
    fn get_neighbor(&self, id: NodeId) -> Result<Option<(NodeId, Edge)>, TreeError>;
    // fn distance(&self, id: NodeId, id2: NodeId) -> Edge;
}

impl TreeHelper for Tree {
    fn get_neighbor(&self, id: NodeId) -> Result<Option<(NodeId, Edge)>, TreeError> {
        let node = self.get(&id)?;
        
        if node.is_root() {
            return Ok(None)
        }

        let parent = node.parent.unwrap();
        let parent_node = self.get(&parent)?;

        let children = parent_node.children.iter().filter(|n| **n != id).collect::<Vec::<&NodeId>>();
        assert!(&children.len() < &2);

        match children.first() {
            Some(child) => {
                let dist = self.get(&parent)?.children.iter().fold(0f64, |acc, child_id| acc + parent_node.get_child_edge(child_id).unwrap());
                Ok(Some((**child, dist)))
            },
            None => Ok(None),
        }
    }

    // fn get_leaves
}

struct LeakageCluster {
    id: NodeId,
    set: HashSet<NodeId>,
    in_: Vec<Leakage>, // Confounding mapping from surrounding 
    out: Vec<Leakage>, // Leakage from this cluster into others
}

impl LeakageCluster {
    pub fn pairwise_distances(&self, tree: &Tree) -> Vec<(NodeId, NodeId, Edge)> {
        let v = Vec::from_iter(&self.set);
        let mut result = Vec::default();

        v.iter().enumerate().for_each(|(idx, source)| {
            let iter = &v[idx+1..].into_iter().map(|target| {
                (**source, **target, tree.get_distance(source, target).unwrap().0.unwrap())
            });
            result.extend(iter.clone());
        });
        
        result
    }
}

fn clean_labels(tree: &mut Tree) {
    let all_nodes = tree.search_nodes(|_| true);
    for nid in &all_nodes {
        
        let node = tree.get_mut(nid).unwrap();
        
        match node.name.as_mut() {
            Some(name) => {
                if name.contains("\"") {
                    let new_label = &name.replace("\"", "");
                    name.clear();
                    name.push_str(&new_label);
                }    
            },
            None => (),
        }
    }
}

fn clean_newick(newick: &str) -> String {
    let single_quotes = newick.chars().filter(|c| *c == '\'').count();
    let double_quotes = newick.chars().filter(|c| *c == '"').count();

    let mut newick_str = String::default();
    if single_quotes > 0 && double_quotes == 0 {
        newick_str = newick.replace("'", "\"");
    }
    newick_str
}

pub fn old_main() {
    println!("Hello, world!");


    let file_path = Path::new("data/trees/bac120_r214.sp_labels.tree");
    let mut newick_str = std::fs::read_to_string(file_path).expect("Cannot read newick-tree from file");

    let single_quotes = newick_str.chars().filter(|c| *c == '\'').count();
    let double_quotes = newick_str.chars().filter(|c| *c == '"').count();

    if single_quotes > 0 && double_quotes == 0 {
        newick_str = newick_str.replace("'", "\"");
        // newick_str = newick_str.replace("'", "");
    }


    let mut tree = match Tree::from_newick(&newick_str) {
        Ok(tree) => tree,
        Err(err) => panic!("{}", err),
    };
    clean_labels(&mut tree);

    eprintln!("Number of leaves: {}", tree.n_leaves());

    
    eprintln!("Number of leaves: {} ... {}", tree.n_leaves(), tree.get_leaf_names().len());

    eprintln!("{:?}", tree.get(&22421));
    eprintln!("{:?}", tree.get(&22430));
    eprintln!("{:?}", tree.get(&22419));

    eprintln!("Common ancestor: {} {} -> {}", 
        22421,
        22430,
        tree.get_common_ancestor(&22421usize, &22430usize).unwrap());


    eprintln!("Number of leaves: {} ... {}", tree.n_leaves(), tree.get_leaf_names().len());
    
    let find = "PUNC01";

    // eprintln!("Id: {:?}", tree.get_node_id(find));
    // eprintln!("Id: {:?}", tree.get_node(find));
    
    let target_nodes = tree.search_nodes(|n| {
        match &n.name {
            Some(name) => name.contains(find),
            None => false,
        }
    });//.iter().map(|id| tree.get(id).unwrap()).collect::<Vec<&Node>>();

    let all_nodes = tree.search_nodes(|_| true);

    let neighbor = tree.get_neighbor(*target_nodes.first().unwrap());
    eprintln!("Neighbor: {:?}", neighbor);

    eprintln!("Nodes {:?}", target_nodes);
    eprintln!("Len all nodes {}", all_nodes.len());


    let mut count = 0;
    for node in &all_nodes {
        let neighbor = tree.get_neighbor(*node);
        // eprintln!("Neighbor: {:?}", neighbor);
        count += 1;
    }
    eprintln!("Count: {}", count);


    // for leaf in &tree.get_leaves() {
    //     let leaf_node = tree.get(leaf).unwrap();
    //     eprintln!("leaf: {:?} -> Neighbor: {:?}", leaf_node.name, tree.get_neighbor(*leaf));
    // }

    let genera = tree.search_nodes(|n| {
        match &n.name {
            Some(name) => name.contains(&"g__"),
            None => false,
        }
    });
    let phyla = tree.search_nodes(|n| {
        n.name.as_ref().map_or(false, |name| name.contains(&"p__"))
    });
    let class = tree.search_nodes(|n| {
        n.name.as_ref().map_or(false, |name| name.contains(&"c__"))
    });
    let order = tree.search_nodes(|n| {
        n.name.as_ref().map_or(false, |name| name.contains(&"o__"))
    });
    let family = tree.search_nodes(|n| {
        n.name.as_ref().map_or(false, |name| name.contains(&"f__"))
    });
    let species = tree.search_nodes(|n| {
        n.name.as_ref().map_or(false, |name| name.contains(&"s__"))
    });




    eprintln!("Phyla len: {}", phyla.len());
    eprintln!("Class len: {}", class.len());
    eprintln!("Order len: {}", order.len());
    eprintln!("Family len: {}", family.len());
    eprintln!("Genera len: {}", genera.len());
    eprintln!("Species len: {}", species.len());
    eprintln!("Leaves len: {}", tree.get_leaves().len());

    for s in &tree.get_leaves() {
        eprintln!("Species name: {:?}", tree.get(s).unwrap());
    }
}

pub fn new_main(newick: String, map: impl AsRef<Path>, leakage_path: impl AsRef<Path>) {
    let newick = clean_newick(&newick);

    let (id2lab, lab2id) = get_labels_map(map);
    
    let mut tree = match Tree::from_newick(&newick) {
        Ok(tree) => tree,
        Err(err) => panic!("{}", err),
    };
    clean_labels(&mut tree);

    /////
    
    let leakage = read_leakage_file(leakage_path);

    eprintln!("Leakage file: {}", leakage.len());

    // Species pair leakage
    let mut species_pair_leakage = HashMap::new();
    for l in &leakage {
        *species_pair_leakage.entry(l.key()).or_insert(0) += 1;
        // eprintln!("From {} To {}", id2lab[l.from], id2lab[l.to] );
    }

    

    let mut sorted_leakage: Vec<(&(usize, usize), &usize)> = species_pair_leakage.iter().collect();
    sorted_leakage.sort_by_key(|e| -(*e.1 as i32));

    sorted_leakage.iter().take(10).for_each(|((t1, t2), events)| {
        eprintln!("{} {} -> {}", id2lab[*t1], id2lab[*t2], events);
    });


    let leakage_counters = get_leakage_counter(&leakage);
    for (id, l) in &leakage_counters {
        eprintln!("{} ({}) -> {}", id2lab[*id], id, l);
    }

}

fn summarize() {
    // Collect command line arguments
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <input_file> <output_file>", args[0]);
        process::exit(1);
    }

    let input_file = &args[1];
    let output_file = &args[2];

    let leakage_summary = read_leakage_counter(input_file);
    
    let mut writer = BufWriter::new(File::create(output_file).unwrap());

    for (id, item) in leakage_summary {
        writer.write_fmt(format_args!("{}\t{}\n", id, item)).expect("Error writing leakage");
    }
}


fn main() {
    // let file_path: &Path = Path::new("data/trees/bac120_r214.sp_labels.tree");
    // let newick_str = std::fs::read_to_string(file_path).expect("Cannot read newick-tree from file");
    // let map_path = "data/maps/genome2tiid.tsv";
    // let leakage_path = "data/leakage_data/61046.bt.summary";

    // new_main(newick_str, map_path, leakage_path);

    summarize();    
}
