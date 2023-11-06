//! Rewrite of back_to_sequence to try improve performance

#![warn(missing_docs)]

/* std use */

/* crate use */

/* project use */

/* mod declaration */
pub mod cli;
pub mod error;
pub mod format;
pub mod io;
pub mod kmer_set;
pub mod tokenizer;

/* pub use */
pub use kmer_set::*;
pub use tokenizer::*;

/// Alias for define Kmer
pub type Kmer = Vec<u8>;
