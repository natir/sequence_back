//! Run a real word benchmark

/* std use */
use std::io::Write as _;

/* crate use */
use anyhow::Context as _;
use clap::Parser as _;

/* mod declaration */
mod common;

/* project use */
use sequence_back::error;

#[derive(clap::Parser, std::fmt::Debug)]
#[clap(
    name = "benches",
    version = "0.1",
    author = "Pierre Marijon <pierre@marijon.fr>"
)]
pub struct Command {
    /// Directory where benchmark data will be write default: benches_data
    #[clap(short = 'd', long = "directory")]
    directory: Option<std::path::PathBuf>,

    /// File where benchmark result will be write default: benches_result.csv
    #[clap(short = 'o', long = "output")]
    output: Option<std::path::PathBuf>,

    /// Generate data
    #[clap(short = 'g', long = "generate", default_value_t = false)]
    generate: bool,

    /// Run sequence_back
    #[clap(short = 's', long = "sequence-back", default_value_t = false)]
    sequence_back: bool,

    /// sequence_back path default: sequence_back
    #[clap(short = 'S', long = "sequence-back-path")]
    sequence_back_path: Option<std::path::PathBuf>,

    /// Run back_to_sequence
    #[clap(short = 'b', long = "back-to-sequence", default_value_t = false)]
    back_to_sequence: bool,

    /// back_to_sequence path default: back_to_sequence
    #[clap(short = 'B', long = "back-to-sequence-path")]
    back_to_sequence_path: Option<std::path::PathBuf>,

    /// Number of run
    #[clap(short = 'n', long = "number-of-run", default_value_t = 10)]
    number_of_run: usize,

    /// Number of reads default: 10000, 100000, 1000000, 10000000, 100000000
    #[clap(short = 'r', long = "reads-numbers", num_args = 1..)]
    reads_numbers: Option<Vec<usize>>,
}

impl Command {
    pub fn directory(&self) -> std::path::PathBuf {
        self.directory
            .clone()
            .unwrap_or(std::path::PathBuf::from("benches_data"))
    }

    pub fn output(&self) -> std::path::PathBuf {
        self.output
            .clone()
            .unwrap_or(std::path::PathBuf::from("benches_result.csv"))
    }

    pub fn generate(&self) -> bool {
        self.generate
    }

    pub fn sequence_back(&self) -> bool {
        self.sequence_back
    }

    pub fn sequence_back_path(&self) -> std::path::PathBuf {
        self.sequence_back_path
            .clone()
            .unwrap_or(std::path::PathBuf::from("sequence_back"))
    }

    pub fn back_to_sequence(&self) -> bool {
        self.back_to_sequence
    }

    pub fn back_to_sequence_path(&self) -> std::path::PathBuf {
        self.back_to_sequence_path
            .clone()
            .unwrap_or(std::path::PathBuf::from("back_to_sequences"))
    }

    pub fn number_of_run(&self) -> usize {
        self.number_of_run
    }

    pub fn reads_numbers(&self) -> Vec<usize> {
        self.reads_numbers
            .clone()
            .unwrap_or(vec![10000, 100000, 1000000, 10000000, 100000000])
    }
}

fn main() -> error::Result<()> {
    // parse cli
    let params = Command::parse();

    stderrlog::new()
        .module(module_path!())
        .verbosity(5)
        .timestamp(stderrlog::Timestamp::Millisecond)
        .init()
        .context("stderrlog already create a logger")?;

    let mut output = std::fs::File::create(params.output())?;
    writeln!(output, "tools,run_index,reads_number,runtime(nanos)",)?;

    if params.generate() {
        std::fs::create_dir_all(params.directory())?;

        let mut rng = common::generator::rng();

        log::info!("Start create reference.");
        let reference = common::generator::fasta(&mut rng, 100000000, 1);
        common::io::write_buffer(&reference, params.directory().join("ref_seq.fasta"))?;
        log::info!("End   create reference.");

        log::info!("Start create kmer.");
        let kmers = common::generator::kmer_from_fasta(&mut rng, &reference, 50, 50000);
        common::io::write_buffer(&kmers, params.directory().join("ref_kmers.fasta"))?;
        log::info!("End   create kmer.");

        for reads_number in params.reads_numbers() {
            log::info!("Start create {} reads.", reads_number);
            let reads =
                common::generator::reads_from_fasta(&mut rng, &reference, 100, reads_number);

            common::io::write_buffer(
                &reads,
                params
                    .directory()
                    .join(format!("reads_{}.fasta", reads_number)),
            )?;
            log::info!("End   create {} reads.", reads_number);
        }
    }

    if params.sequence_back() {
        for reads_number in params.reads_numbers() {
            for index in 0..params.number_of_run() {
                log::info!("Start run sequence_back on {}", reads_number);
                let mut cmd = assert_cmd::Command::new(params.sequence_back_path());

                cmd.args([
                    "--input-kmers",
                    params
                        .directory()
                        .join("ref_kmers.fasta")
                        .as_os_str()
                        .to_str()
                        .unwrap(),
                    "--input-sequences",
                    params
                        .directory()
                        .join(format!("reads_{}.fasta", reads_number))
                        .as_os_str()
                        .to_str()
                        .unwrap(),
                    "--output-sequences",
                    params
                        .directory()
                        .join(format!("filtered_{}.fasta", reads_number))
                        .as_os_str()
                        .to_str()
                        .unwrap(),
                    "-k",
                    "31",
                    "--output-kmers",
                    params
                        .directory()
                        .join(format!("counted_kmer_{}.csv", reads_number))
                        .as_os_str()
                        .to_str()
                        .unwrap(),
                ]);

                let start = std::time::Instant::now();
                cmd.assert().success();
                let runtime = std::time::Instant::now().duration_since(start);

                writeln!(
                    output,
                    "sequence_back,{},{},{}",
                    index,
                    reads_number,
                    runtime.as_nanos()
                )?;
                log::info!("End   run sequence_back on {}", reads_number);
            }
        }
    }

    if params.back_to_sequence() {
        for reads_number in params.reads_numbers() {
            for index in 0..params.number_of_run() {
                log::info!("Start run back_to_sequence on {}", reads_number);
                let mut cmd = assert_cmd::Command::new(params.back_to_sequence_path());

                cmd.args([
                    "--in-kmers",
                    params
                        .directory()
                        .join("ref_kmers.fasta")
                        .as_os_str()
                        .to_str()
                        .unwrap(),
                    "--in-sequences",
                    params
                        .directory()
                        .join(format!("reads_{}.fasta", reads_number))
                        .as_os_str()
                        .to_str()
                        .unwrap(),
                    "--out-sequences",
                    params
                        .directory()
                        .join(format!("filtered_{}.fasta", reads_number))
                        .as_os_str()
                        .to_str()
                        .unwrap(),
                    "-k",
                    "31",
                    "--out-kmers",
                    params
                        .directory()
                        .join(format!("counted_kmer_{}.csv", reads_number))
                        .as_os_str()
                        .to_str()
                        .unwrap(),
                ]);

                let start = std::time::Instant::now();
                cmd.assert().success();
                let runtime = std::time::Instant::now().duration_since(start);

                writeln!(
                    output,
                    "back_to_sequence,{},{},{}",
                    index,
                    reads_number,
                    runtime.as_nanos()
                )?;
                log::info!("End   run back_to_sequence on {}", reads_number);
            }
        }
    }

    Ok(())
}
