//! Error struct of project sequence_back

/* crate use */
use anyhow;
use thiserror;

/// Enum to manage error
#[derive(std::fmt::Debug, thiserror::Error)]
pub enum Error {
    /// Input stream isn't fasta or fastq format
    #[error("Input stream isn't in fasta or fastq format")]
    NotFastaOrFastq,

    /// Set option query-reverse useless without stranded option
    #[error("Set option query-reverse useless without stranded option")]
    QueryReverseWithoutStranded,

    /// Error in logging system configuration
    #[error(transparent)]
    Log(#[from] log::SetLoggerError),
}

/// Alias of result
pub type Result<T> = anyhow::Result<T>;
