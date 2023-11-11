//! Benchmark forward kmer selection

/* std use */

/* crate use */
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

/* mod declaration */
mod common;

/* project use */
use common::*;
use sequence_back::{KmerCounter, KmerCounterTrait};

fn number_of_reads(c: &mut Criterion) {
    let mut g = c.benchmark_group("number_of_reads");

    g.sample_size(50);
    g.warm_up_time(std::time::Duration::from_secs(1));

    let mut rng = generator::rng();

    for pow in 10..17 {
        let number_of_sequence = 2u64.pow(pow);

        let fasta = common::generator::fasta(&mut rng, 50, number_of_sequence);

        g.bench_with_input(
            BenchmarkId::from_parameter(number_of_sequence),
            &fasta,
            |b, fasta| {
                b.iter(|| {
                    black_box(KmerCounter::from_fasta_stream(
                        fasta.as_slice(),
                        31,
                        false,
                        false,
                        8192,
                    ))
                })
            },
        );
    }
}

fn kmer_size(c: &mut Criterion) {
    let mut g = c.benchmark_group("kmer_size");

    g.sample_size(40);
    g.warm_up_time(std::time::Duration::from_secs(1));

    let mut rng = generator::rng();

    for kmer_size in (10..97).step_by(2) {
        let fasta = common::generator::fasta(&mut rng, 100, 8192);

        g.bench_with_input(
            BenchmarkId::from_parameter(kmer_size),
            &fasta,
            |b, fasta| {
                b.iter(|| {
                    black_box(KmerCounter::from_fasta_stream(
                        fasta.as_slice(),
                        kmer_size,
                        false,
                        false,
                        8192,
                    ))
                })
            },
        );
    }
}

fn buffer_size(c: &mut Criterion) {
    let mut g = c.benchmark_group("buffer_size");

    g.sample_size(40);
    g.warm_up_time(std::time::Duration::from_secs(1));

    let mut rng = generator::rng();

    for pow in 10..17 {
        let buffer_size = 2u64.pow(pow);
        let fasta = common::generator::fasta(&mut rng, 50, 65536);

        g.bench_with_input(
            BenchmarkId::from_parameter(buffer_size),
            &fasta,
            |b, fasta| {
                b.iter(|| {
                    black_box(KmerCounter::from_fasta_stream(
                        fasta.as_slice(),
                        31,
                        false,
                        false,
                        buffer_size,
                    ))
                })
            },
        );
    }
}

fn number_of_thread(c: &mut Criterion) {
    let mut g = c.benchmark_group("number_of_thread");

    g.sample_size(40);
    g.warm_up_time(std::time::Duration::from_secs(1));

    let mut rng = generator::rng();

    for number_of_threads in (2..17).step_by(2) {
        let fasta = common::generator::fasta(&mut rng, 50, 65536);

        g.bench_with_input(
            BenchmarkId::from_parameter(number_of_threads),
            &fasta,
            |b, fasta| {
                b.iter(|| {
                    std::env::set_var("RAYON_NUM_THREADS", format!("{}", number_of_threads));

                    black_box(KmerCounter::from_fasta_stream(
                        fasta.as_slice(),
                        31,
                        false,
                        false,
                        8192,
                    ))
                })
            },
        );
    }
}

fn kmer_counter_init(c: &mut Criterion) {
    number_of_reads(c);
    kmer_size(c);
    buffer_size(c);
    number_of_thread(c);
}

criterion_group!(benches, kmer_counter_init);

criterion_main!(benches);
