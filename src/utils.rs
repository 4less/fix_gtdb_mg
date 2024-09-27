use std::{fs::File, io::{BufRead, BufReader, Read, Seek, SeekFrom}, path::Path};

use flate2::bufread::GzDecoder;

/// Returns an iterator over the lines of a given file, handling both plain text and gzipped files.
pub fn file_lines<P: AsRef<Path>>(path: P) -> std::io::Result<Box<dyn Iterator<Item = std::io::Result<String>>>> {
    let file = File::open(&path)?;
    // Check if the file is gzipped by looking at its first two bytes
    let mut buf_reader = BufReader::new(file);
    let mut first_two_bytes = [0u8; 2];
    
    buf_reader.read_exact(&mut first_two_bytes)?;

    // Reset the reader to start again
    buf_reader.seek(SeekFrom::Start(0));
    let reader: Box<dyn BufRead> = if first_two_bytes == [0x1F, 0x8B] {
        // It's a gzipped file
        let decoder = GzDecoder::new(buf_reader);
        Box::new(BufReader::new(decoder))
    } else {
        // It's a regular text file
        Box::new(buf_reader)
    };
    
    Ok(Box::new(reader.lines()))
}
