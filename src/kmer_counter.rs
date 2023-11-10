//! Specialize an HashSet to get intresting kmer.

/* std use */
use std::io::Read as _;

/* crate use */
use ahash;
use atomic_counter::AtomicCounter;
use rayon::prelude::*;

/* project use */
use crate::error;
use crate::format;
use crate::get_tokenizer;
use crate::io;
use crate::Kmer;

macro_rules! from_stream {
    ($name:ident, $reader:ty, $populate:ident, $doc:ident) => {
        /// Build an KmerCounter from a $doc stream
        fn $name<R>(
            input: R,
            kmer_size: u64,
            stranded: bool,
            reverse: bool,
            record_buffer_len: u64,
        ) -> error::Result<KmerCounter>
        where
            R: std::io::BufRead,
        {
            let hasher = ahash::random_state::RandomState::with_seed(42);
            let mut hash_set: KmerCounter = KmerCounter::with_capacity_and_hasher(65536, hasher);
            let mut reader = <$reader>::new(input);
            let mut iter = reader.records();
            let mut records = Vec::with_capacity(record_buffer_len as usize);

            let mut end = true;
            while end {
                end = io::$populate(&mut iter, &mut records, record_buffer_len);

                hash_set.par_extend(
                    records
                        .par_iter()
                        .map(|record| {
                            let norm_seq = record.sequence().as_ref().to_ascii_uppercase();
                            get_tokenizer(&norm_seq, kmer_size, stranded, reverse)
                                .unwrap()
                                .map(|x| (x, atomic_counter::RelaxedCounter::default()))
                                .collect::<Vec<(Vec<u8>, atomic_counter::RelaxedCounter)>>()
                        })
                        .flatten(),
                );
            }

            Ok(hash_set)
        }
    };
}

/// Trait to add functionalty to AHashSet
pub trait KmerCounterTrait
where
    Self: Default,
{
    /// Build an KmerCounter from stream
    fn from_stream<R>(
        mut input: R,
        kmer_size: u64,
        stranded: bool,
        reverse: bool,
        record_buffer_len: u64,
    ) -> error::Result<KmerCounter>
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

    from_stream!(
        from_fasta_stream,
        noodles::fasta::Reader<R>,
        fasta_populate_buffer,
        fasta
    );

    from_stream!(
        from_fastq_stream,
        noodles::fastq::Reader<R>,
        fastq_populate_buffer,
        fastq
    );

    /// Write
    fn to_csv<W>(kmerset: KmerCounter, output: &mut W) -> error::Result<()>
    where
        W: std::io::Write,
    {
        writeln!(output, "kmer,count")?;
        for (key, value) in kmerset {
            writeln!(
                output,
                "{},{}",
                std::str::from_utf8(&key).unwrap(),
                value.get()
            )?;
        }

        Ok(())
    }
}

/// Alias for AHashSet<Kmer>
pub type KmerCounter = ahash::AHashMap<Kmer, atomic_counter::RelaxedCounter>;

impl KmerCounterTrait for KmerCounter {}

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
        let kmerset = KmerCounter::from_stream(FASTA_FILE, 4, true, false, 2)?;

        let mut kmers = kmerset.keys().cloned().collect::<Vec<Vec<u8>>>();
        let mut truth = KMERS.iter().map(|a| a.to_vec()).collect::<Vec<Vec<u8>>>();

        kmers.sort();
        truth.sort();

        assert_eq!(kmers, truth);

        let kmerset = KmerCounter::from_stream(FASTQ_FILE, 4, true, false, 2)?;

        let mut kmers = kmerset.keys().cloned().collect::<Vec<Vec<u8>>>();
        let mut truth = KMERS.iter().map(|a| a.to_vec()).collect::<Vec<Vec<u8>>>();

        kmers.sort();
        truth.sort();

        assert_eq!(kmers, truth);

        Ok(())
    }

    #[test]
    fn from_fasta_stream() -> error::Result<()> {
        let kmerset = KmerCounter::from_fasta_stream(FASTA_FILE, 4, true, false, 2)?;

        let mut kmers = kmerset.keys().cloned().collect::<Vec<Vec<u8>>>();
        let mut truth = KMERS.iter().map(|a| a.to_vec()).collect::<Vec<Vec<u8>>>();

        kmers.sort();
        truth.sort();

        assert_eq!(kmers, truth);

        Ok(())
    }

    #[test]
    fn from_fastq_stream() -> error::Result<()> {
        let kmerset = KmerCounter::from_fastq_stream(FASTQ_FILE, 4, true, false, 2)?;

        let mut kmers = kmerset.keys().cloned().collect::<Vec<Vec<u8>>>();
        let mut truth = KMERS.iter().map(|a| a.to_vec()).collect::<Vec<Vec<u8>>>();

        kmers.sort();
        truth.sort();

        assert_eq!(kmers, truth);

        Ok(())
    }

    #[test]
    fn to_csv() -> error::Result<()> {
        let kmerset = KmerCounter::from_fastq_stream(FASTQ_FILE, 4, true, false, 2)?;
        let mut output = Vec::new();

        KmerCounter::to_csv(kmerset, &mut output)?;

        let mut result = output
            .split(|c| *c == b'\n')
            .map(|x| x.to_vec())
            .collect::<Vec<Vec<u8>>>();
        result.sort();

        assert_eq!(
            result,
            vec![
                b"".to_vec(),
                b"AAAA,0".to_vec(),
                b"ACTG,0".to_vec(),
                b"GGGG,0".to_vec(),
                b"GTCA,0".to_vec(),
                b"kmer,count".to_vec(),
            ]
        );

        Ok(())
    }
}
