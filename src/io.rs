//! Specialize an HashSet to get intresting kmer.

/* std use */
use std::io::Read as _;

/* crate use */
use atomic_counter::AtomicCounter;
use rayon::prelude::*;

/* project use */
use crate::error;
use crate::format;
use crate::get_tokenizer;
use crate::kmer_counter;

macro_rules! populate_buffer {
    ($name:ident, $iterator_ty:ty, $records_ty:ty) => {
        /// Populate record buffer with content of iterator
        pub(crate) fn $name<R>(
            iter: &mut $iterator_ty,
            records: &mut $records_ty,
            record_buffer: u64,
        ) -> bool
        where
            R: std::io::BufRead,
        {
            records.clear();

            for i in 0..record_buffer {
                if let Some(Ok(record)) = iter.next() {
                    records.push(record);
                } else {
                    records.truncate(i as usize);
                    return false;
                }
            }

            true
        }
    };
}

populate_buffer!(
    fasta_populate_buffer,
    noodles::fasta::reader::Records<'_, R>,
    Vec<noodles::fasta::Record>
);

populate_buffer!(
    fastq_populate_buffer,
    noodles::fastq::reader::Records<'_, R>,
    Vec<noodles::fastq::Record>
);

#[allow(clippy::too_many_arguments)]
/// Filter reads by kmer
pub fn filter_reads<R, W>(
    mut input: R,
    output: W,
    kmer_counter: &mut kmer_counter::KmerCounter,
    kmer_size: u64,
    min_threshold: f64,
    max_threshold: f64,
    stranded: bool,
    reverse: bool,
    record_buffer_len: u64,
) -> error::Result<()>
where
    R: std::io::BufRead,
    W: std::io::Write,
{
    match format::ReadsFormat::detect(&mut input)? {
        format::ReadsFormat::Fasta => filter_fasta_reads(
            std::io::Cursor::new([b'>']).chain(input),
            output,
            kmer_counter,
            kmer_size,
            min_threshold,
            max_threshold,
            stranded,
            reverse,
            record_buffer_len,
        ),
        format::ReadsFormat::Fastq => filter_fastq_reads(
            std::io::Cursor::new([b'@']).chain(input),
            output,
            kmer_counter,
            kmer_size,
            min_threshold,
            max_threshold,
            stranded,
            reverse,
            record_buffer_len,
        ),
    }
}

macro_rules! filter_reads {
    ( $type:ident, $name:ident, $reader_type:ty, $writer_type:ty, $populate:ident, $record_ty:ty) => {
        #[allow(clippy::too_many_arguments)]
        /// Filter $type reads by kmer
        pub(crate) fn $name<R, W>(
            input: R,
            output: W,
            kmer_counter: &mut kmer_counter::KmerCounter,
            kmer_size: u64,
            min_threshold: f64,
            max_threshold: f64,
            stranded: bool,
            reverse: bool,
            record_buffer_len: u64,
        ) -> error::Result<()>
        where
            R: std::io::BufRead,
            W: std::io::Write,
        {
            let mut reader = <$reader_type>::new(input);
            let mut writer = <$writer_type>::new(output);
            let mut iter = reader.records();
            let mut records = Vec::with_capacity(record_buffer_len as usize);

            let mut end = true;
            while end {
                log::info!("Start read record");
                end = $populate(&mut iter, &mut records, record_buffer_len);
                log::info!("End read record");
                log::info!("Populate buffer with {} record", records.len());

                let filterd_records: $record_ty = records
                    .par_iter()
                    .filter_map(|record| {
                        let norm_seq = record.sequence().as_ref().to_ascii_uppercase();

                        let valid_kmer = get_tokenizer(&norm_seq, kmer_size, stranded, reverse)
                            .unwrap()
                            .filter(|kmer| kmer_counter.contains_key(kmer))
                            .inspect(|kmer| {
                                kmer_counter.get(kmer).unwrap().inc();
                            })
                            .count() as f64;

                        let total_kmer: f64 =
                            record.sequence().len() as f64 - kmer_size as f64 + 1.0;
                        let ratio: f64 = valid_kmer / total_kmer * 100.0;

                        if ratio > min_threshold && ratio < max_threshold {
                            Some(record)
                        } else {
                            None
                        }
                    })
                    .collect();
                log::info!("{} filtred record", filterd_records.len());

                log::info!("Start write record");
                for record in filterd_records {
                    writer.write_record(record)?;
                }
                log::info!("End write record");
            }

            Ok(())
        }
    };
}

filter_reads!(
    fasta,
    filter_fasta_reads,
    noodles::fasta::Reader<R>,
    noodles::fasta::Writer<W>,
    fasta_populate_buffer,
    Vec<&noodles::fasta::Record>
);

filter_reads!(
    fastq,
    filter_fastq_reads,
    noodles::fastq::Reader<R>,
    noodles::fastq::Writer<W>,
    fastq_populate_buffer,
    Vec<&noodles::fastq::Record>
);

#[cfg(test)]
mod tests {
    use super::*;

    const FASTA: &[u8] = b">random_seq 0
AGAAGGCCCCGACTTTCGGGGTGGGAAGTCTTTTGTTCTCACGTTACACCGCCTGCGCAATAGTAGGGCGAAAAATCGCC
>random_seq 1
AAFGGCTGCTCTCCAGCTAATACCATCTGAGCTTTGGCTCAGGTGAAGGCACGACACCAGAAATCAACGGCGAGGCGGTT
>random_seq 2
CATCAGCTACGTGGCGGCCTGAGCAACGGCTCGCTTTGAGTAGTCAATGGATCTGGGTATGGCAGTCATTTTCCCGCTTT
>random_seq 3
CGGCGCTTACTCTATGTGAGCGACGTGGAGAAAGATCCGGAATGACCTAGATAGATTTGCAAGCTTGGACAACTACGGAG
>random_seq 4
GAGCGCAGAAGACCACCAGAAGGATCTCAAGAAGGACCAGGTCAGAGTTAGGCATATTCCCGTTCGAGGCTGATATTGGC
";
    const FASTA_KMER: [&[u8; 15]; 3] = [b"AGCTAATACCATCTG", b"GTCAGAGTTAGGCAT", b"ACCAGAAATCAACGG"];

    mod fasta {
        use super::*;

        #[test]
        fn populate_buffer() {
            let file = std::io::Cursor::new(FASTA);

            let mut reader = noodles::fasta::Reader::new(file);
            let mut iter = reader.records();
            let mut records = Vec::with_capacity(3);

            assert!(fasta_populate_buffer(&mut iter, &mut records, 3));

            assert_eq!(
                records,
                vec![
                    noodles::fasta::Record::new(
                        noodles::fasta::record::Definition::new("random_seq", Some("0".to_string())),
                        noodles::fasta::record::Sequence::from(b"AGAAGGCCCCGACTTTCGGGGTGGGAAGTCTTTTGTTCTCACGTTACACCGCCTGCGCAATAGTAGGGCGAAAAATCGCC".to_vec())
                    ),
                    noodles::fasta::Record::new(
                        noodles::fasta::record::Definition::new("random_seq", Some("1".to_string())),
                        noodles::fasta::record::Sequence::from(b"AAFGGCTGCTCTCCAGCTAATACCATCTGAGCTTTGGCTCAGGTGAAGGCACGACACCAGAAATCAACGGCGAGGCGGTT".to_vec())
                    ),
                    noodles::fasta::Record::new(
                        noodles::fasta::record::Definition::new("random_seq", Some("2".to_string())),
                        noodles::fasta::record::Sequence::from(b"CATCAGCTACGTGGCGGCCTGAGCAACGGCTCGCTTTGAGTAGTCAATGGATCTGGGTATGGCAGTCATTTTCCCGCTTT".to_vec())
                    ),
                ]
            );

            assert!(!fasta_populate_buffer(&mut iter, &mut records, 3));

            assert_eq!(
                records,
		vec![
                    noodles::fasta::Record::new(
                        noodles::fasta::record::Definition::new("random_seq", Some("3".to_string())),
                        noodles::fasta::record::Sequence::from(b"CGGCGCTTACTCTATGTGAGCGACGTGGAGAAAGATCCGGAATGACCTAGATAGATTTGCAAGCTTGGACAACTACGGAG".to_vec())
                    ),
                    noodles::fasta::Record::new(
                        noodles::fasta::record::Definition::new("random_seq", Some("4".to_string())),
                        noodles::fasta::record::Sequence::from(b"GAGCGCAGAAGACCACCAGAAGGATCTCAAGAAGGACCAGGTCAGAGTTAGGCATATTCCCGTTCGAGGCTGATATTGGC".to_vec())
                ),]
            );

            assert!(!fasta_populate_buffer(&mut iter, &mut records, 3));

            assert_eq!(records, vec![]);
        }

        #[test]
        fn filter_reads() -> error::Result<()> {
            let mut output = Vec::new();
            let writer = Box::new(std::io::BufWriter::new(std::io::Cursor::new(&mut output)));
            let mut counter = FASTA_KMER
                .iter()
                .map(|x| (x.to_vec(), atomic_counter::RelaxedCounter::default()))
                .collect();

            filter_fasta_reads(
                FASTA,
                writer,
                &mut counter,
                15,
                0.0,
                100.0,
                true,
                false,
                8192,
            )?;

            assert_eq!(
                output,
                b">random_seq 1
AAFGGCTGCTCTCCAGCTAATACCATCTGAGCTTTGGCTCAGGTGAAGGCACGACACCAGAAATCAACGGCGAGGCGGTT
>random_seq 4
GAGCGCAGAAGACCACCAGAAGGATCTCAAGAAGGACCAGGTCAGAGTTAGGCATATTCCCGTTCGAGGCTGATATTGGC
"
            );

            let mut output = Vec::new();
            let writer = Box::new(std::io::BufWriter::new(std::io::Cursor::new(&mut output)));
            let mut counter = FASTA_KMER
                .iter()
                .map(|x| (x.to_vec(), atomic_counter::RelaxedCounter::default()))
                .collect();

            filter_fasta_reads(
                FASTA,
                writer,
                &mut counter,
                15,
                2.0,
                100.0,
                true,
                false,
                8192,
            )?;

            assert_eq!(
                output,
                b">random_seq 1
AAFGGCTGCTCTCCAGCTAATACCATCTGAGCTTTGGCTCAGGTGAAGGCACGACACCAGAAATCAACGGCGAGGCGGTT
"
            );

            Ok(())
        }
    }

    const FASTQ: &[u8] = b"@random_seq 0
CCATACGAAAGTTCGAGGGACAAGGATGATGCCTCGGATCTTATGAAGTATCCCGTTAACCAGTTATAGGCTAATGTAGG
+
XmaM3IM.kloVw92/SAHr>#2!6oMnn=32cNa/8L58{9P~CS2$'HRrs8Tt:{j=c)O<oqX+]weMd3_PK%|{
@random_seq 1
TAAGCTGACTGGTCGTAAAGTGAGCTCGCATCCATGCTTCGCATCAAGTGAGATAAGAGTTACCACTGTGGCGTTCGTGA
+
'@Y]WtoA%)u307)ZMC_G+spyv`$pCDv?7q4Nsg)um;Zf7GqsOY?TR%^{lsJlC@;bZQ)Z$okSDcTZ$p(N
@random_seq 2
GCATTAATACTATCCGACGCCCGGTGGCTCTTTCCAAACGTTGCAATGTATATAGCGTACACCTGGAGTTAAGACAGTCG
+
9QP<77+dj<@%nI&,7HW/^~J/yV7}s4WPFG@S|neKu.7So'cuhI^{NwA5<Oh7G4|<?D@CSRQTmTG`I3)J
@random_seq 3
TCAGGGGTGCCTCATAGGGTCCCCTCCCCCGGGCCCGTTGTTATAACTTTTGCACGTAACCCTATGGAACTGCTAAGCGT
+
&DQU-A7!.`:v*+bJKb0oGJ&<(}LpGCBYRKq4O}}2w(rffZ)-x<c})1YQb)@M_zvozOo(|@M1%aBIg%Y{
@random_seq 4
GCACGCAGAACAAAGGCGCGACTCCCCAAGGCTACAACTTTATGGCCGCTCCGTCGAGTGGGCCATACAGGATATCGTAA
+
RSO8jKE0{jQ<_^AS`}Vp$5/<a]G}%qx;Eian<+vqt.WHXoKSV]9ASb%$|E*tO#>l$I6bx_/o.h&_r%Z~
";

    const FASTQ_KMER: [&[u8; 15]; 3] = [b"CGAAAGTTCGAGGGA", b"TGCCTCATAGGGTCC", b"GCCCGTTGTTATAAC"];

    mod fastq {
        use super::*;

        #[test]
        fn populate_buffer() {
            let file = std::io::Cursor::new(FASTQ);

            let mut reader = noodles::fastq::Reader::new(file);
            let mut iter = reader.records();
            let mut records = Vec::with_capacity(3);

            assert!(fastq_populate_buffer(&mut iter, &mut records, 3));

            assert_eq!(
                records,
                vec![
                    noodles::fastq::Record::new(
                        noodles::fastq::record::Definition::new("random_seq", "0"),
                        "CCATACGAAAGTTCGAGGGACAAGGATGATGCCTCGGATCTTATGAAGTATCCCGTTAACCAGTTATAGGCTAATGTAGG",
                        "XmaM3IM.kloVw92/SAHr>#2!6oMnn=32cNa/8L58{9P~CS2$'HRrs8Tt:{j=c)O<oqX+]weMd3_PK%|{",
                    ),
                    noodles::fastq::Record::new(
                        noodles::fastq::record::Definition::new("random_seq", "1"),
                        "TAAGCTGACTGGTCGTAAAGTGAGCTCGCATCCATGCTTCGCATCAAGTGAGATAAGAGTTACCACTGTGGCGTTCGTGA",
                        "'@Y]WtoA%)u307)ZMC_G+spyv`$pCDv?7q4Nsg)um;Zf7GqsOY?TR%^{lsJlC@;bZQ)Z$okSDcTZ$p(N",
                    ),
                    noodles::fastq::Record::new(
                        noodles::fastq::record::Definition::new("random_seq", "2"),
                        "GCATTAATACTATCCGACGCCCGGTGGCTCTTTCCAAACGTTGCAATGTATATAGCGTACACCTGGAGTTAAGACAGTCG",
                        "9QP<77+dj<@%nI&,7HW/^~J/yV7}s4WPFG@S|neKu.7So'cuhI^{NwA5<Oh7G4|<?D@CSRQTmTG`I3)J",
                    ),
                ]
            );

            assert!(!fastq_populate_buffer(&mut iter, &mut records, 3));

            assert_eq!(
                records,
                vec![
		    noodles::fastq::Record::new(
                        noodles::fastq::record::Definition::new("random_seq", "3"),
                        "TCAGGGGTGCCTCATAGGGTCCCCTCCCCCGGGCCCGTTGTTATAACTTTTGCACGTAACCCTATGGAACTGCTAAGCGT",
                        "&DQU-A7!.`:v*+bJKb0oGJ&<(}LpGCBYRKq4O}}2w(rffZ)-x<c})1YQb)@M_zvozOo(|@M1%aBIg%Y{",
                    ),
                    noodles::fastq::Record::new(
                        noodles::fastq::record::Definition::new("random_seq", "4"),
                        "GCACGCAGAACAAAGGCGCGACTCCCCAAGGCTACAACTTTATGGCCGCTCCGTCGAGTGGGCCATACAGGATATCGTAA",
                        "RSO8jKE0{jQ<_^AS`}Vp$5/<a]G}%qx;Eian<+vqt.WHXoKSV]9ASb%$|E*tO#>l$I6bx_/o.h&_r%Z~",
                    ),
                ]
            );

            assert!(!fastq_populate_buffer(&mut iter, &mut records, 3));

            assert_eq!(records, vec![]);
        }

        #[test]
        fn filter_reads() -> error::Result<()> {
            let mut output = Vec::new();
            let writer = Box::new(std::io::BufWriter::new(std::io::Cursor::new(&mut output)));
            let mut counter = FASTQ_KMER
                .iter()
                .map(|x| (x.to_vec(), atomic_counter::RelaxedCounter::default()))
                .collect();

            filter_fastq_reads(
                FASTQ,
                writer,
                &mut counter,
                15,
                0.0,
                100.0,
                true,
                false,
                8192,
            )?;

            assert_eq!(
                output,
                b"@random_seq 0
CCATACGAAAGTTCGAGGGACAAGGATGATGCCTCGGATCTTATGAAGTATCCCGTTAACCAGTTATAGGCTAATGTAGG
+
XmaM3IM.kloVw92/SAHr>#2!6oMnn=32cNa/8L58{9P~CS2$'HRrs8Tt:{j=c)O<oqX+]weMd3_PK%|{
@random_seq 3
TCAGGGGTGCCTCATAGGGTCCCCTCCCCCGGGCCCGTTGTTATAACTTTTGCACGTAACCCTATGGAACTGCTAAGCGT
+
&DQU-A7!.`:v*+bJKb0oGJ&<(}LpGCBYRKq4O}}2w(rffZ)-x<c})1YQb)@M_zvozOo(|@M1%aBIg%Y{
"
            );

            let mut output = Vec::new();
            let writer = Box::new(std::io::BufWriter::new(std::io::Cursor::new(&mut output)));
            let mut counter = FASTQ_KMER
                .iter()
                .map(|x| (x.to_vec(), atomic_counter::RelaxedCounter::default()))
                .collect();

            filter_fastq_reads(
                FASTQ,
                writer,
                &mut counter,
                15,
                3.0,
                100.0,
                true,
                false,
                8192,
            )?;

            assert_eq!(
                output,
                b"@random_seq 3
TCAGGGGTGCCTCATAGGGTCCCCTCCCCCGGGCCCGTTGTTATAACTTTTGCACGTAACCCTATGGAACTGCTAAGCGT
+
&DQU-A7!.`:v*+bJKb0oGJ&<(}LpGCBYRKq4O}}2w(rffZ)-x<c})1YQb)@M_zvozOo(|@M1%aBIg%Y{
"
            );

            Ok(())
        }
    }

    #[test]
    fn filter_reads_() -> error::Result<()> {
        let mut output = Vec::new();
        let writer = Box::new(std::io::BufWriter::new(std::io::Cursor::new(&mut output)));
        let mut counter = FASTA_KMER
            .iter()
            .map(|x| (x.to_vec(), atomic_counter::RelaxedCounter::default()))
            .collect();

        filter_reads(
            FASTA,
            writer,
            &mut counter,
            15,
            0.0,
            100.0,
            true,
            false,
            8192,
        )?;

        assert_eq!(
            output,
            b">random_seq 1
AAFGGCTGCTCTCCAGCTAATACCATCTGAGCTTTGGCTCAGGTGAAGGCACGACACCAGAAATCAACGGCGAGGCGGTT
>random_seq 4
GAGCGCAGAAGACCACCAGAAGGATCTCAAGAAGGACCAGGTCAGAGTTAGGCATATTCCCGTTCGAGGCTGATATTGGC
"
        );

        let mut output = Vec::new();
        let writer = Box::new(std::io::BufWriter::new(std::io::Cursor::new(&mut output)));
        let mut counter = FASTQ_KMER
            .iter()
            .map(|x| (x.to_vec(), atomic_counter::RelaxedCounter::default()))
            .collect();

        filter_reads(
            FASTQ,
            writer,
            &mut counter,
            15,
            0.0,
            100.0,
            true,
            false,
            8192,
        )?;

        assert_eq!(
            output,
            b"@random_seq 0
CCATACGAAAGTTCGAGGGACAAGGATGATGCCTCGGATCTTATGAAGTATCCCGTTAACCAGTTATAGGCTAATGTAGG
+
XmaM3IM.kloVw92/SAHr>#2!6oMnn=32cNa/8L58{9P~CS2$'HRrs8Tt:{j=c)O<oqX+]weMd3_PK%|{
@random_seq 3
TCAGGGGTGCCTCATAGGGTCCCCTCCCCCGGGCCCGTTGTTATAACTTTTGCACGTAACCCTATGGAACTGCTAAGCGT
+
&DQU-A7!.`:v*+bJKb0oGJ&<(}LpGCBYRKq4O}}2w(rffZ)-x<c})1YQb)@M_zvozOo(|@M1%aBIg%Y{
"
        );

        Ok(())
    }
}
