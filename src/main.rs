//! Rewrite of sequence_back to try improve performance

#![warn(missing_docs)]

/* std use */

/* crate use */
use anyhow::Context as _;
use clap::Parser as _;

/* project use */
use sequence_back::cli;
use sequence_back::error;
use sequence_back::io;
use sequence_back::{KmerCounter, KmerCounterTrait};

fn main() -> error::Result<()> {
    // parse cli
    let params = cli::Command::parse();

    // Setup logger
    stderrlog::new()
        .module(module_path!())
        .quiet(params.quiet())
        .verbosity(params.verbosity())
        .timestamp(params.timestamp())
        .init()
        .context("stderrlog already create a logger")?;

    log::info!("Start build hash_set");
    let mut kmer_set = KmerCounter::from_stream(
        params.input_kmers()?,
        params.kmer_size(),
        params.stranded(),
        params.query_reverse(),
        params.record_buffer(),
    )?;
    log::info!("End build hash_set");

    log::info!("Start filter reads");
    io::filter_reads(
        params.input_sequences()?,
        params.output_sequences()?,
        &mut kmer_set,
        params.kmer_size(),
        params.min_threshold(),
        params.max_threshold(),
        params.stranded(),
        params.query_reverse(),
        params.record_buffer(),
    )?;
    log::info!("End filter reads");

    log::info!("Start write kmer counts");
    KmerCounter::to_csv(kmer_set, &mut params.output_kmers()?)?;
    log::info!("End write kmer counts");

    Ok(())
}
