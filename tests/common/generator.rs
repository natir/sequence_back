//! Define function to generate sequence value use for test and benchmarking

/* std use */

/* crate use */
use rand::prelude::*;
use rayon::slice::ParallelSliceMut;

/* project use */
use crate::common::constant;
use sequence_back::kmer_counter::{KmerCounter, KmerCounterTrait};

/// Generate a RNG with constant::SEED
pub fn rng() -> rand::rngs::StdRng {
    rand::rngs::StdRng::from_seed(constant::SEED)
}

/// Generate a random DNA sequence with constant::SEQUENCE_ALPHABET
pub fn seq(rng: &mut rand::rngs::StdRng, seq_length: u64) -> Vec<u8> {
    (0..seq_length)
        .map(|_| *constant::SEQUENCE_ALPHABET.choose(rng).unwrap())
        .collect::<Vec<u8>>()
}

/// Generate a random DNA sequence with constant::QUALITY_ALPHABET
pub fn quality(rng: &mut rand::rngs::StdRng, seq_length: u64) -> Vec<u8> {
    (0..seq_length)
        .map(|_| *constant::QUALITY_ALPHABET.choose(rng).unwrap())
        .collect::<Vec<u8>>()
}

/// Generate a random in ram fasta with RNG
pub fn fasta(rng: &mut rand::rngs::StdRng, seq_length: u64, seq_number: u64) -> Vec<u8> {
    let mut output = Vec::with_capacity(
        (
            seq_length * seq_number // sequence space
		+ seq_number * 3 // '>' and jump line space
		+ (seq_number.checked_ilog10().unwrap_or(0) as u64 + 1) * seq_number
            // sequence id space
        ) as usize,
    );

    for index in 0..seq_number {
        // Header
        output.extend(b">");
        output.extend(index.to_string().as_bytes());
        output.extend(b"\n");
        // Sequence
        output.extend(seq(rng, seq_length));
        output.extend(b"\n");
    }

    output
}

/// Generate a random in ram fastq with RNG
pub fn fastq(rng: &mut rand::rngs::StdRng, seq_length: u64, seq_number: u64) -> Vec<u8> {
    let mut output = Vec::with_capacity(
        (
            seq_length * seq_number // sequence space
            + seq_length * seq_number // quality space
            + seq_number * 6 // '@' '+' and jump line space
	+ (seq_number.checked_ilog10().unwrap_or(0) as u64 + 1) * seq_number
            // sequence id space
        ) as usize,
    );

    for index in 0..seq_number {
        // Header
        output.extend(b"@");
        output.extend(index.to_string().as_bytes());
        output.extend(b"\n");
        // Sequence
        output.extend(seq(rng, seq_length));
        output.extend(b"\n");
        // Plus
        output.extend(b"+\n");
        // Quality
        output.extend(quality(rng, seq_length));
        output.extend(b"\n");
    }

    output
}

#[allow(dead_code)]
/// Generate kmer from fasta buffer
pub fn kmer_from_fasta(
    rng: &mut rand::rngs::StdRng,
    in_buffer: &[u8],
    kmers_size: u64,
    kmers_number: u64,
) -> Vec<u8> {
    let mut kmerset =
        KmerCounter::from_fasta_stream(in_buffer, kmers_size, false, false, 8192).unwrap();

    let mut out_buffer = Vec::with_capacity(
        (kmers_size * kmers_number // sequence space
            + kmers_number * 3 // '>' and jump line space
	    + (kmers_number.checked_ilog10().unwrap_or(0) as u64 + 1) * kmers_number) // sequence id space
     as usize,
    );

    let mut kmers: Vec<Vec<u8>> = kmerset.drain().map(|(kmer, _count)| kmer).collect();
    kmers.par_sort_unstable();

    for (index, kmer) in kmers
        .choose_multiple(rng, kmers_number as usize)
        .enumerate()
    {
        // Header
        out_buffer.extend(b">");
        out_buffer.extend(index.to_string().as_bytes());
        out_buffer.extend(b"\n");
        // Sequence
        out_buffer.extend(kmer);
        out_buffer.extend(b"\n");
    }

    out_buffer
}

#[allow(dead_code)]
/// Generate kmer from fastq buffer
pub fn kmer_from_fastq(
    rng: &mut rand::rngs::StdRng,
    in_buffer: &[u8],
    kmers_size: u64,
    kmers_number: u64,
) -> Vec<u8> {
    let mut kmerset =
        KmerCounter::from_fastq_stream(in_buffer, kmers_size, false, false, 8192).unwrap();

    let mut out_buffer = Vec::with_capacity(
        (kmers_size * kmers_number // sequence space
            + kmers_size * kmers_number // quality space
            + kmers_number * 6 // '@' '+' and jump line space
	    + (kmers_number.checked_ilog10().unwrap_or(0) as u64 + 1) * kmers_number) // sequence id space
     as usize,
    );

    let mut kmers: Vec<Vec<u8>> = kmerset.drain().map(|(kmer, _count)| kmer).collect();
    kmers.par_sort_unstable();

    for (index, kmer) in kmers
        .choose_multiple(rng, kmers_number as usize)
        .enumerate()
    {
        // Header
        out_buffer.extend(b">");
        out_buffer.extend(index.to_string().as_bytes());
        out_buffer.extend(b"\n");
        // Sequence
        out_buffer.extend(kmer);
        out_buffer.extend(b"\n");
    }

    out_buffer
}

#[allow(dead_code)]
/// Generate reads from fasta buffer
pub fn reads_from_fasta(
    rng: &mut rand::rngs::StdRng,
    in_buffer: &[u8],
    reads_size: usize,
    reads_number: usize,
) -> Vec<u8> {
    let mut reader = noodles::fasta::Reader::new(in_buffer);

    let mut out_buffer = Vec::with_capacity(
        reads_size * reads_number // sequence space
            + reads_number * 3 // '>' and jump line space
	    + (reads_number.checked_ilog10().unwrap_or(0) as usize + 1) * reads_number, // sequence id space
    );

    let sequences = reader
        .records()
        .map(|record| record.unwrap().sequence().as_ref().to_vec())
        .collect::<Vec<Vec<u8>>>();

    for index in 0..reads_number {
        let seq = sequences.choose(rng).unwrap();
        let pos = rng.gen_range(0..(seq.len() - reads_size));

        // Header
        out_buffer.extend(b">");
        out_buffer.extend(index.to_string().as_bytes());
        out_buffer.extend(b"\n");

        // Sequence
        out_buffer.extend(&seq[pos..(pos + reads_size)]);
        out_buffer.extend(b"\n");
    }

    out_buffer
}

#[allow(dead_code)]
/// Generate reads from fasta buffer
pub fn reads_from_fastq(
    rng: &mut rand::rngs::StdRng,
    in_buffer: &[u8],
    reads_size: usize,
    reads_number: usize,
) -> Vec<u8> {
    let mut reader = noodles::fasta::Reader::new(in_buffer);

    let mut out_buffer = Vec::with_capacity(
        reads_size * reads_number // sequence space
            + reads_size * reads_number // quality space
            + reads_number * 6 // '@' '+' and jump line space
	    + (reads_number.checked_ilog10().unwrap_or(0) as usize + 1) * reads_number, // sequence id space
    );

    let sequences = reader
        .records()
        .map(|record| record.unwrap().sequence().as_ref().to_vec())
        .collect::<Vec<Vec<u8>>>();

    for index in 0..reads_number {
        let seq = sequences.choose(rng).unwrap();
        let pos = rng.gen_range(0..(seq.len() - reads_size));

        // Header
        out_buffer.extend(b"@");
        out_buffer.extend(index.to_string().as_bytes());
        out_buffer.extend(b"\n");
        // Sequence
        out_buffer.extend(&seq[pos..(pos + reads_size)]);
        out_buffer.extend(b"\n");
        // Plus
        out_buffer.extend(b"+\n");
        // Quality
        out_buffer.extend(quality(rng, reads_size as u64));
        out_buffer.extend(b"\n");
    }

    out_buffer
}
