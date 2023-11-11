<h1 style="text-align: center;">sequence_back</h1>

![Test](https://github.com/natir/sequence_back/workflows/Test/badge.svg)
![Lints](https://github.com/natir/sequence_back/workflows/Lints/badge.svg)
![MSRV](https://github.com/natir/sequence_back/workflows/MSRV/badge.svg)
[![CodeCov](https://codecov.io/gh/natir/sequence_back/branch/main/graph/badge.svg)](https://codecov.io/gh/natir/sequence_back)
[![Documentation](https://github.com/natir/sequence_back/workflows/Documentation/badge.svg)](https://natir.github.io/sequence_back/sequence_back)
[![License](https://img.shields.io/badge/license-MIT-green)](https://github.com/natir/sequence_back/blob/master/LICENSE)


Rewrite of [back_to_sequence](https://github.com/pierrepeterlongo/back_to_sequences) to try improve performance

## Installation

### With cargo

```bash
cargo install --git https://github.com/natir/sequence_back.git
```

### From source

```bash
git clone https://github.com/natir/sequence_back.git
cd sequence_back
cargo install --path .
```

## Usage

### Help

```
Rewrite of back_to_sequence to try improve performance

Usage: sequence_back [OPTIONS] --input-kmers <INPUT_KMERS> --output-kmers <OUTPUT_KMERS>

Options:
  -i, --input-sequences <INPUT_SEQUENCES>
          Input fasta or fastq file containing the original sequences
  -I, --input-kmers <INPUT_KMERS>
          Input fasta or fastq file containing the original kmers
  -o, --output-sequences <OUTPUT_SEQUENCES>
          Output file containing the filtred sequences
  -O, --output-kmers <OUTPUT_KMERS>
          Output file containing the filtred kmers
  -k, --kmer-size <KMER_SIZE>
          Size of the kmers to index and search [default: 31]
  -m, --min-threshold <MIN_THRESHOLD>
          Output sequences are those whose ratio of indexed kmers is in ]min_threshold; max_threshold] Minimal threshold of the ratio  (%) of kmers that must be found in a sequence to keep it (default 0%). Thus by default, if no kmer is found in a sequence, it is not output [default: 0]
  -M, --max-threshold <MAX_THRESHOLD>
          Output sequences are those whose ratio of indexed kmers is in ]min_threshold; max_threshold] Maximal threshold of the ratio (%) of kmers that must be found in a sequence to keep it (default 100%). Thus by default, there is no limitation on the maximal number of kmers found in a sequence [default: 100]
  -s, --stranded
          Used original kmer strand (else canonical kmers are considered)
  -r, --query-reverse
          Query the reverse complement of reads. Useless without the --stranded option
  -b, --record-buffer <RECORD_BUFFER>
          Control size of record buffer [default: 8192]
  -t, --threads <THREADS>
          Number of theard use 0 use all avaible core, default value 0 [default: 0]
  -q, --quiet
          Silence all output
  -v, --verbosity...
          Verbose mode (-v, -vv, -vvv, etc)
  -T, --timestamp <TS>
          Timestamp (sec, ms, ns, none)
  -h, --help
          Print help
  -V, --version
          Print version
```

### Examples

#### Basic

```bash
sequences_back --input-kmers compacted_kmers.fasta --input-sequences reads.fasta --output-sequences filtered_reads.fasta  --output-kmers counted_kmers.csv
```

The `filtered_reads.fasta` file contains the original sequences (here reads) from reads.fasta that contain at least one of the kmers in `compacted_kmers.fasta`.
The file `counted_kmers.csv` contains for each kmer in `compacted_kmers.fasta` the number of times it was found in `filtered_reads.fasta`.

#### Using filters

```bash
sequences_back --input-kmers compacted_kmers.fasta --input-sequences reads.fasta --output-sequences filtered_reads.fasta  --output-kmers counted_kmers.csv --min-threshold 50 --max-threshold 70
```

In this case only sequeces from `reads.fasta` that have more than 50% and at most 70% of their kmers in `compacted_kmers.fasta` are output.

## Benchmark

```
# Generate data
cargo run --example benches --release -- -g
# Run sequence_back
cargo run --example benches --release -- -s
# Run sequence_back and back_to_sequence
cargo run --example benches --release -- -s -b
```

**Warning**: by default this benches assume `sequence_back` and `back_to_sequence` are install in your `PATH`, you can use `--sequence-back-path` and `--back-to-sequence-path` to set path of these tools.

## Minimum supported Rust version

Currently the minimum supported Rust version is 1.73.
