use bitm::{BitAccess, BitVec};
use bitnuc_mismatch::generate_mismatches;
use clap::Parser;
use serde_json::json;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::vec::Vec;
use zip::ZipArchive;

#[derive(Parser, Debug)]
#[command(name = "sequence-matcher")]
#[command(about = "Matches TSV sequences against a reference set from a zip file")]
struct Args {
    /// Path to the zip file containing known permit list sequences
    #[arg(short, long)]
    zip_file: PathBuf,

    /// Path to the TSV file to search
    #[arg(short, long)]
    tsv_file: PathBuf,
}

fn read_sequences_from_zip(
    zip_path: &PathBuf,
) -> Result<(Box<[u64]>, usize), Box<dyn std::error::Error>> {
    let file = File::open(zip_path)?;
    let mut archive = ZipArchive::new(file)?;
    let mut parent_sequences = Vec::with_capacity(3_000_000);
    let mut sequences = Box::<[u64]>::with_zeroed_bits(4_usize.pow(16_u32));

    // Read the first file in the archive (adjust if needed)
    let mut zip_file = archive.by_index(0)?;
    let reader = BufReader::new(&mut zip_file);

    for line in reader.lines() {
        let sequence = line?;
        let trimmed = sequence.trim();
        if trimmed.len() == 16 {
            parent_sequences.push(trimmed.to_string());
        }
    }

    let parent_scalars: Vec<u64> = parent_sequences
        .into_iter()
        .map(|seq| bitnuc::as_2bit(seq.as_bytes()).unwrap())
        .collect();

    let mut mm_buffer = Vec::with_capacity(3_usize.pow(16_u32));

    let mut nset = 0;
    for p in &parent_scalars {
        if !sequences.get_bit(*p as usize) {
            nset += 1;
            sequences.set_bit(*p as usize);
        }
        generate_mismatches(*p, 16, &mut mm_buffer)?;

        for mm in &mm_buffer {
            if !sequences.get_bit(*mm as usize) {
                nset += 1;
                sequences.set_bit(*mm as usize);
            }
        }
    }

    Ok((sequences, nset))
}

fn process_tsv(
    tsv_path: &PathBuf,
    reference_set: &Box<[u64]>,
) -> Result<(usize, usize), Box<dyn std::error::Error>> {
    let file = File::open(tsv_path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    let mut matched_count = 0;
    let mut total_count = 0;

    while let Some(line) = lines.next() {
        let l = line?;

        // Read sequence line
        if let Some(bc_entry) = l.split_ascii_whitespace().skip(1).next() {
            total_count += 1;

            // Check if umi is in 1 mismatch set
            let bc = unsafe { bitnuc::as_2bit(bc_entry.get_unchecked(0..16).as_bytes())? };

            if reference_set.get_bit(bc as usize) {
                matched_count += 1;
            }
        }
    }

    Ok((total_count, matched_count))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    let (reference_set, _nset) = read_sequences_from_zip(&args.zip_file)?;

    eprintln!("Processing TSV file...");
    let (tot, matched) = process_tsv(&args.tsv_file, &reference_set)?;
    let matched_frac = (matched as f64) / (tot as f64);

    let output = json!({
        "total_barcodes": tot,
        "matched_barcodes": matched,
        "matched_frac": matched_frac
    });

    println!("{}", serde_json::to_string_pretty(&output)?);
    Ok(())
}
