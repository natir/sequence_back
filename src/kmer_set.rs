//! Specialize an HashSet to get intresting kmer.

/* std use */
use std::io::Read as _;

/* crate use */
use ahash;
use rayon::prelude::*;

/* project use */
use crate::error;
use crate::format;
use crate::get_tokenizer;
use crate::io;
use crate::Kmer;

/// Trait to add functionalty to AHashSet
pub trait KmerSetTrait
where
    Self: Default,
{
    /// Build an KmerSet from stream
    fn from_stream<R>(
        mut input: R,
        kmer_size: u64,
        stranded: bool,
        reverse: bool,
        record_buffer_len: u64,
    ) -> error::Result<KmerSet>
    where
        R: std::io::BufRead,
    {
        match format::ReadsFormat::detect(&mut input)? {
            format::ReadsFormat::Fasta => Self::from_fasta_stream(
                std::io::Cursor::new([b'>']).chain(input),
                kmer_size,
                stranded,
                reverse,
                record_buffer_len,
            ),
            format::ReadsFormat::Fastq => Self::from_fastq_stream(
                std::io::Cursor::new([b'@']).chain(input),
                kmer_size,
                stranded,
                reverse,
                record_buffer_len,
            ),
        }
    }

    /// Build an KmerSet from a fasta stream
    fn from_fasta_stream<R>(
        input: R,
        kmer_size: u64,
        stranded: bool,
        reverse: bool,
        record_buffer_len: u64,
    ) -> error::Result<KmerSet>
    where
        R: std::io::BufRead,
    {
        let hasher = ahash::random_state::RandomState::with_seed(42);
        let mut hash_set: KmerSet = KmerSet::with_capacity_and_hasher(65536, hasher);
        let mut reader = noodles::fasta::Reader::new(input);
        let mut iter = reader.records();
        let mut records = Vec::with_capacity(record_buffer_len as usize);

        let mut end = true;
        while end {
            end = io::fasta_populate_buffer(&mut iter, &mut records, record_buffer_len);

            hash_set.par_extend(
                records
                    .par_iter()
                    .map(|record| {
                        let norm_seq = record.sequence().as_ref().to_ascii_uppercase();
                        get_tokenizer(&norm_seq, kmer_size, stranded, reverse)
                            .unwrap()
                            .collect::<Vec<Vec<u8>>>()
                    })
                    .flatten(),
            );
        }

        Ok(hash_set)
    }

    /// Build an KmerSet from a fastq stream
    fn from_fastq_stream<R>(
        input: R,
        kmer_size: u64,
        stranded: bool,
        reverse: bool,
        record_buffer_len: u64,
    ) -> error::Result<KmerSet>
    where
        R: std::io::BufRead,
    {
        let hasher = ahash::random_state::RandomState::with_seed(42);
        let mut hash_set: KmerSet = KmerSet::with_capacity_and_hasher(65536, hasher);
        let mut reader = noodles::fastq::Reader::new(input);
        let mut iter = reader.records();
        let mut records = Vec::with_capacity(record_buffer_len as usize);

        let mut end = true;
        while end {
            end = io::fastq_populate_buffer(&mut iter, &mut records, record_buffer_len);

            hash_set.par_extend(
                records
                    .par_iter()
                    .map(|record| {
                        let norm_seq = record.sequence().as_ref().to_ascii_uppercase();
                        get_tokenizer(&norm_seq, kmer_size, stranded, reverse)
                            .unwrap()
                            .collect::<Vec<Vec<u8>>>()
                    })
                    .flatten(),
            );
        }

        Ok(hash_set)
    }
}

/// Alias for AHashSet<Kmer>
pub type KmerSet = ahash::AHashSet<Kmer>;

impl KmerSetTrait for KmerSet {}

#[cfg(test)]
mod tests {
    use super::*;

    const FASTA_FILE: &[u8] = b">1
ACTG
>2
GTCA
>3
AAAA
>4
GGGG";

    const FASTQ_FILE: &[u8] = b"@1
ACTG
+
!!!!
@2
GTCA
+
!!!!
@3
AAAA
+
!!!!
@4
GGGG
+
!!!!";

    const KMERS: &[&[u8]] = &[b"ACTG", b"GTCA", b"AAAA", b"GGGG"];

    #[test]
    fn from_stream() -> error::Result<()> {
        let kmerset = KmerSet::from_stream(FASTA_FILE, 4, true, false, 2)?;

        assert_eq!(kmerset, KMERS.iter().map(|a| a.to_vec()).collect());

        let kmerset = KmerSet::from_stream(FASTQ_FILE, 4, true, false, 2)?;

        assert_eq!(kmerset, KMERS.iter().map(|a| a.to_vec()).collect());

        Ok(())
    }

    #[test]
    fn from_fasta_stream() -> error::Result<()> {
        let kmerset = KmerSet::from_fasta_stream(FASTA_FILE, 4, true, false, 2)?;

        assert_eq!(kmerset, KMERS.iter().map(|a| a.to_vec()).collect());

        Ok(())
    }

    #[test]
    fn from_fastq_stream() -> error::Result<()> {
        let kmerset = KmerSet::from_fastq_stream(FASTQ_FILE, 4, true, false, 2)?;

        assert_eq!(kmerset, KMERS.iter().map(|a| a.to_vec()).collect());

        Ok(())
    }
}
