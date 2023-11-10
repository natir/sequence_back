//! Functional test on sequence_back

/* std use */

/* crate use */

/* mod declaration */
mod common;

/* project use */
use common::*;

#[test]
fn help() -> Result<(), anyhow::Error> {
    let mut cmd = assert_cmd::Command::cargo_bin("sequence_back").unwrap();

    cmd.args(["-h"]);

    let truth: &[u8] = b"Rewrite of back_to_sequence to try improve performance

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
";

    let assert = cmd.assert();

    assert.success().stdout(truth);

    Ok(())
}

#[test]
fn default_fasta() -> Result<(), anyhow::Error> {
    let mut cmd = assert_cmd::Command::cargo_bin("sequence_back").unwrap();
    let mut rng = generator::rng();

    let temp_dir = tempfile::tempdir()?;
    let temp_path = temp_dir.path();
    let kmers_in_path = temp_path.join("kmers_in.fasta");
    let kmers_out_path = temp_path.join("kmers_out.fasta");

    let reads = generator::fasta(&mut rng, 150, 100);
    let kmers = generator::kmer_from_fasta(&mut rng, &reads, 31, 5);

    io::write_buffer(&kmers, &kmers_in_path)?;

    let truth: &[u8] = b">9
ctGgccgtGAgAGacGtGCAGTAGaGagttGcatCGaaGGtaTGtGagCGaCAAGtgtgtTtcAaACttGACGttgatAA
ACtgattatctCTCcattCgatCAaGaTatcaTTagtATTcgtaGgTaCgAtcGaAGcGCacGaattTaA
>24
TGCcTgagcacccactgtTTtAaCgttcAttttaCccCaGCCGGTcTgtTCACccGACAaGgTTTCGTGGGtgcTCcAag
attgcaTAGaAGGGtgAtGgaGCCgAaGCGCgCaGaATgTaCcACcccTCAaggGcAaatatGgaTGGcT
>38
cgGCCCTtTaTgCAcTcgGAATtCCGAtcTTcCtaCCGTgGCgCcAAgtgcCttACGacAGaGTcAgaCcATAtCgttAA
AGcGCtAGtGGtAggtcaGtGAaACTtCgcTTccTTtCCTGACtGTCTTgGcCgTggCAaTttCtCatAa
>51
AcGtAGGGGatgTCgggaCaGTATATtccTCATaTTccgCGTCaaCGaCTatgtgACTtcTcGtTccCAcTTACTAacCC
ATCACtgATtCggggTCctCGGCtgTggatgtCccgtacAgGctaCGAgTTgGtCtgCTggAaCCtCtCG
>92
TGCAgcccaTCCAcCgaGGaTATAtcCcTaTAcGCGGgGcTtGCCtGacTTCgaCcAAtCTTggTTAgCgtATatAaAcG
aAATTactTaaGtgCctTAtgttcaccaCCCAgaTgaccTGATctatTATgcCCAtTgcAgAtcGTCAAT
";

    cmd.args([
        "-I",
        kmers_in_path.as_os_str().to_str().unwrap(),
        "-O",
        kmers_out_path.as_os_str().to_str().unwrap(),
    ])
    .write_stdin(reads);

    let assert = cmd.assert();

    assert.success().stdout(truth);

    Ok(())
}

#[test]
fn default_fastq() -> Result<(), anyhow::Error> {
    let mut cmd = assert_cmd::Command::cargo_bin("sequence_back").unwrap();
    let mut rng = generator::rng();

    let temp_dir = tempfile::tempdir()?;
    let temp_path = temp_dir.path();
    let kmers_in_path = temp_path.join("kmers_in.fasta");
    let kmers_out_path = temp_path.join("kmers_out.fasta");

    let reads = generator::fastq(&mut rng, 150, 100);

    let kmers = generator::kmer_from_fastq(&mut rng, &reads, 31, 5);

    io::write_buffer(&kmers, &kmers_in_path)?;

    let truth: &[u8] = b"@42
cCACcCgAtaaGTCactAcTCTGtaAcCGAAGTtatAcATAGTaACTcgATTaCTtctgcttCagCCtaCttTtaGTGGtATGcAGgCaAgCcTGTTggGCtCCTGggCcGTcAAACTttcTCgAACGgAaAtcggTAcGtcGCaTTgcT
+
570960571469771674334657910663679884647962396416299990599594236190526900279695565765851925046961179917632156482025062988276584859798522178769489435585
@59
cTTCgaTgCcaaGcgGtAaggaaGtgTacAatggTaggATtacCTcgCCATAtACggAgCaCtttTgAaTcaCgAtcTGCgaatgTtGaaGGccAtCtTGAccCCAAtGCcgcTGgctAaTGTtGTCTcgAcATGgaAacCaTGTcaCTA
+
628709127791110608280812496630402970244531206686010665946303076313437084189268307552413020627659035990044177172449417344211650649881276332383424356714
@76
cCCAtctaCAGcATCCtagAGAAAaacGaaGgGagaGCaAaAatAcTtcgTttgAggGATaTttgTGTAaATCTCCggCTcgaTcgGctcTAGccaCcGGTTCgtCaAGgtTggTACGcCtCctAcATGaAACtttgTTtAttGaTCtaT
+
480153694100453102292790226804147064578716149458468250542976903642209518506573763063528934282936658737281813965900325604160274286667043319091459160303
@86
GACtGAGGtCggaccgcgacGtCgacTgctTAaGcGaAgggaTCggAcaACctaaAtTctCgatgtcAtAAGccAagaACcAgACGtGAGCcgcGaCgcGATtccgCaCgCAGGGggTcataGattTACCcAaCgACCTcACaCACCaAg
+
455058777669714361639995673590285097835256620977948381331811150165433005043318167533573129615624547217962736149939168738435264143055522171752934182368
";

    cmd.args([
        "-I",
        kmers_in_path.as_os_str().to_str().unwrap(),
        "-O",
        kmers_out_path.as_os_str().to_str().unwrap(),
    ])
    .write_stdin(reads);

    let assert = cmd.assert();

    assert.success().stdout(truth);

    Ok(())
}

#[test]
fn k15() -> Result<(), anyhow::Error> {
    let mut cmd = assert_cmd::Command::cargo_bin("sequence_back").unwrap();
    let mut rng = generator::rng();

    let temp_dir = tempfile::tempdir()?;
    let temp_path = temp_dir.path();
    let kmers_in_path = temp_path.join("kmers_in.fasta");
    let kmers_out_path = temp_path.join("kmers_out.fasta");

    let reads = generator::fasta(&mut rng, 150, 100);
    let kmers = generator::kmer_from_fasta(&mut rng, &reads, 15, 5);

    io::write_buffer(&kmers, &kmers_in_path)?;

    let truth: &[u8] = b">43
CTGtGGTTATgAgGgaaTCaAgTcaGtAtCTATTcgGaACTATgGGCtAgccATCcCggtCCAgCaaCAccctAtgcggt
aCAtGgaggCCgCTaTCGTtCaCGTgTaaGAggaGggtCTCaCaCCtctAcgGCATaGAtAaAgTgTctc
>46
CAcTAtacTGGAtCaggggatcCttgCACtaCaTAggTTctCtAAGAtTATACtCcCCagTataCtGTCGtATccAtcTg
aAGGCcaCgAAcCgaaTTTcTgatgaGACtccgCaAGGtaaAAgaaAGGATCTatgaatGcgGtgTCgAC
>69
acaCAgGaTtagtTGTAcaTGTCctTTcCaTctttttaAatTacTacTtggTAacacattaaGatGTtCctaTTccgAgc
tGACTaGActgtaGcgCccgcTGGCAGAaCcGAtaCCACcGCACTCaCcTgAaTCcCAaggctgTatcgC
>91
tTaCTACtcCACcACctAtaAtTcGctAccTgTTTGgCGCtCGGcCTCgCaCgtCcAaGaGCtaGacatgGAGGTtAAaG
GCTctACtTTttccAACGaGcATCggaGTcttCcGcctAcTacgtcCactgGgGCCTGCaaatggtccAc
";

    cmd.args([
        "-I",
        kmers_in_path.as_os_str().to_str().unwrap(),
        "-O",
        kmers_out_path.as_os_str().to_str().unwrap(),
        "-k",
        "15",
    ])
    .write_stdin(reads);

    let assert = cmd.assert();

    assert.success().stdout(truth);

    Ok(())
}
