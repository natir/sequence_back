//! Command Line Interface declaration of project sequence_back

/* std use */

/* crate use */

/* project use */
use crate::error::Result;

/// Rewrite of back_to_sequence to try improve performance
#[derive(clap::Parser, std::fmt::Debug)]
#[clap(
    name = "sequence_back",
    version = "0.1",
    author = "Pierre Marijon <pierre@marijon.fr>"
)]
pub struct Command {
    /* Specifique option */
    /// Input fasta or fastq file containing the original sequences
    #[clap(short = 'i', long = "input-sequences")]
    input_sequences: Option<std::path::PathBuf>,

    /// Input fasta or fastq file containing the original kmers
    #[clap(short = 'I', long = "input-kmers")]
    input_kmers: std::path::PathBuf,

    /// Output file containing the filtred sequences.
    #[clap(short = 'o', long = "output-sequences")]
    output_sequences: Option<std::path::PathBuf>,

    /// Output file containing the filtred kmers.
    #[clap(short = 'O', long = "output-kmers")]
    output_kmers: std::path::PathBuf,

    /// Size of the kmers to index and search
    #[clap(short = 'k', long = "kmer-size", default_value_t = 31)]
    kmer_size: u64,

    /// Output sequences are those whose ratio of indexed kmers is in ]min_threshold; max_threshold]
    /// Minimal threshold of the ratio  (%) of kmers that must be found in a sequence to keep it (default 0%).
    /// Thus by default, if no kmer is found in a sequence, it is not output.
    #[clap(short = 'm', long = "min-threshold", default_value_t = 0.0)]
    min_threshold: f64,

    /// Output sequences are those whose ratio of indexed kmers is in ]min_threshold; max_threshold]
    /// Maximal threshold of the ratio (%) of kmers that must be found in a sequence to keep it (default 100%).
    /// Thus by default, there is no limitation on the maximal number of kmers found in a sequence.
    #[clap(short = 'M', long = "max-threshold", default_value_t = 100.0)]
    max_threshold: f64,

    /// Used original kmer strand (else canonical kmers are considered)
    #[arg(short = 's', long = "stranded", default_value_t = false)]
    stranded: bool,

    /// Query the reverse complement of reads. Useless without the --stranded option
    #[arg(short = 'r', long = "query-reverse", default_value_t = false)]
    query_reverse: bool,

    /// Control size of record buffer
    #[arg(short = 'b', long = "record-buffer", default_value_t = 8192)]
    record_buffer: u64,

    /* General option */
    /// Silence all output
    #[clap(short = 'q', long = "quiet")]
    quiet: bool,

    /// Verbose mode (-v, -vv, -vvv, etc)
    #[clap(short = 'v', long = "verbosity", action = clap::ArgAction::Count)]
    verbosity: u8,

    /// Timestamp (sec, ms, ns, none)
    #[clap(short = 'T', long = "timestamp")]
    ts: Option<stderrlog::Timestamp>,
}

impl Command {
    /// Get sequences reader
    pub fn input_sequences(&self) -> Result<Box<std::io::BufReader<dyn std::io::Read>>> {
        if let Some(path) = &self.input_sequences {
            Ok(Box::new(std::io::BufReader::new(
                niffler::get_reader(Box::new(std::fs::File::open(path)?))?.0,
            )))
        } else {
            Ok(Box::new(std::io::BufReader::new(std::io::stdin())))
        }
    }

    /// Get kmers reader
    pub fn input_kmers(&self) -> Result<Box<std::io::BufReader<dyn std::io::Read>>> {
        Ok(Box::new(std::io::BufReader::new(
            niffler::get_reader(Box::new(std::fs::File::open(&self.input_kmers)?))?.0,
        )))
    }

    /// Get sequences writer
    pub fn output_sequences(&self) -> Result<Box<std::io::BufWriter<dyn std::io::Write>>> {
        if let Some(path) = &self.output_sequences {
            Ok(Box::new(std::io::BufWriter::new(Box::new(
                std::fs::File::create(path)?,
            ))))
        } else {
            Ok(Box::new(std::io::BufWriter::new(std::io::stdout())))
        }
    }

    /// Get kmers writer
    pub fn output_kmers(&self) -> Result<Box<std::io::BufReader<dyn std::io::Read>>> {
        Ok(Box::new(std::io::BufReader::new(std::fs::File::create(
            &self.output_kmers,
        )?)))
    }

    /// Get kmer size
    pub fn kmer_size(&self) -> u64 {
        self.kmer_size
    }

    /// Get minimal threshold
    pub fn min_threshold(&self) -> f64 {
        self.min_threshold
    }

    /// Get maximal threshold
    pub fn max_threshold(&self) -> f64 {
        self.max_threshold
    }

    /// Get stranded
    pub fn stranded(&self) -> bool {
        self.stranded
    }

    /// Get record buffer length
    pub fn record_buffer(&self) -> u64 {
        self.record_buffer
    }

    /// Get query_reverse
    pub fn query_reverse(&self) -> bool {
        self.query_reverse
    }

    /// Get verbosity level
    pub fn verbosity(&self) -> usize {
        self.verbosity as usize
    }

    /// Get quiet
    pub fn quiet(&self) -> bool {
        self.quiet
    }

    /// Get timestamp granularity
    pub fn timestamp(&self) -> stderrlog::Timestamp {
        self.ts.unwrap_or(stderrlog::Timestamp::Off)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use clap::Parser as _;
    use std::io::Write as _;

    #[test]
    fn default() -> Result<()> {
        let tmp_path = tempfile::tempdir()?;
        let mut in_kmer = std::fs::File::create(tmp_path.path().join("in_kmers.fasta"))?;
        let _out_kmer = std::fs::File::create(tmp_path.path().join("out_kmers.fasta"))?;

        in_kmer.write_all(b"ACTGA")?;

        let params = Command::parse_from([
            "sequence_back".as_ref(),
            "-I".as_ref(),
            tmp_path.path().join("in_kmers.fasta").as_os_str(),
            "-O".as_ref(),
            tmp_path.path().join("out_kmers.fasta").as_os_str(),
        ]);

        assert!(params.input_sequences().is_ok());
        assert!(params.input_kmers().is_ok());
        assert!(params.output_sequences().is_ok());
        assert!(params.output_kmers().is_ok());

        assert_eq!(params.kmer_size(), 31);
        assert_eq!(params.min_threshold(), 0.0);
        assert_eq!(params.max_threshold(), 100.0);
        assert!(!params.stranded());
        assert!(!params.query_reverse());
        assert_eq!(params.record_buffer, 8192);

        assert!(!params.quiet());
        assert_eq!(params.verbosity(), 0);
        assert!(matches!(params.timestamp(), stderrlog::Timestamp::Off));

        Ok(())
    }

    #[test]
    fn not_default() -> Result<()> {
        let tmp_path = tempfile::tempdir()?;
        let mut in_seq = std::fs::File::create(tmp_path.path().join("in_sequences.fasta"))?;
        let mut in_kmer = std::fs::File::create(tmp_path.path().join("in_kmers.fasta"))?;
        let _out_seq = std::fs::File::create(tmp_path.path().join("out_sequences.fasta"))?;
        let _out_kmer = std::fs::File::create(tmp_path.path().join("out_kmers.fasta"))?;

        in_seq.write_all(b"ACTGA")?;
        in_kmer.write_all(b"ACTGA")?;

        let params = Command::parse_from([
            "sequence_back".as_ref(),
            "-i".as_ref(),
            tmp_path.path().join("in_sequences.fasta").as_os_str(),
            "-I".as_ref(),
            tmp_path.path().join("in_kmers.fasta").as_os_str(),
            "-o".as_ref(),
            tmp_path.path().join("out_sequences.fasta").as_os_str(),
            "-O".as_ref(),
            tmp_path.path().join("out_kmers.fasta").as_os_str(),
            "-k".as_ref(),
            "63".as_ref(),
            "-m".as_ref(),
            "20.0".as_ref(),
            "-M".as_ref(),
            "80.0".as_ref(),
            "-s".as_ref(),
            "-r".as_ref(),
            "-b".as_ref(),
            "42".as_ref(),
            "-q".as_ref(),
            "-vvvv".as_ref(),
            "-T".as_ref(),
            "ns".as_ref(),
        ]);

        assert!(params.input_sequences().is_ok());
        assert!(params.input_kmers().is_ok());
        assert!(params.output_sequences().is_ok());
        assert!(params.output_kmers().is_ok());

        assert_eq!(params.kmer_size(), 63);
        assert_eq!(params.min_threshold(), 20.0);
        assert_eq!(params.max_threshold(), 80.0);
        assert!(params.stranded());
        assert!(params.query_reverse());
        assert_eq!(params.record_buffer(), 42);

        assert!(params.quiet());
        assert_eq!(params.verbosity(), 4);
        assert!(matches!(
            params.timestamp(),
            stderrlog::Timestamp::Nanosecond
        ));

        Ok(())
    }
}
