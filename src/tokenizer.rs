//! Define iterator over kmer of sequences

/* std use */

/* crate use */

/* project use */
use crate::error;

const fn comp_lookup_table() -> [u8; 256] {
    const fn complement(nuc: u8) -> u8 {
        match nuc {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b'T',
            b't' => b'A',
            b'c' => b'G',
            b'g' => b'C',
            b'n' => b'N',
            n => n,
        }
    }

    let mut i = 0;
    let mut lookup = [0; 256];
    while i != 256 {
        lookup[i] = complement(i as u8);
        i += 1;
    }

    lookup
}

const COMP: [u8; 256] = comp_lookup_table();

/// Create the good tokenizer
pub fn get_tokenizer<'a>(
    seq: &'a [u8],
    kmer_size: u64,
    stranded: bool,
    reverse: bool,
) -> error::Result<Box<dyn Iterator<Item = Vec<u8>> + 'a>> {
    match (stranded, reverse) {
        (false, false) => Ok(Box::new(Canonical::new(seq, kmer_size))),
        (true, false) => Ok(Box::new(Forward::new(seq, kmer_size))),
        (true, true) => Ok(Box::new(Reverse::new(seq, kmer_size))),
        (false, true) => Err(error::Error::QueryReverseWithoutStranded.into()),
    }
}

/// Struct that perform tokenization in forward mode
struct Forward<'a> {
    seq: &'a [u8],
    kmer_size: usize,
    start: usize,
}

impl<'a> Forward<'a> {
    /// Create a new forward tokenizer
    pub fn new(seq: &'a [u8], kmer_size: u64) -> Self {
        Self {
            seq,
            kmer_size: kmer_size as usize,
            start: 0,
        }
    }
}

impl<'a> Iterator for Forward<'a> {
    type Item = Vec<u8>;

    /// Create a new kmer
    fn next(&mut self) -> Option<Self::Item> {
        if self.start + self.kmer_size > self.seq.len() {
            None
        } else {
            let kmer = self.seq[self.start..self.start + self.kmer_size].to_vec();
            self.start += 1;
            Some(kmer)
        }
    }
}

/// Struct that perform tokenization in reverse mode
struct Reverse {
    seq: Vec<u8>,
    kmer_size: usize,
    start: usize,
}

impl Reverse {
    /// Create a new reverse tokenizer
    pub fn new(seq: &[u8], kmer_size: u64) -> Self {
        let revcomp = seq
            .iter()
            .rev()
            .map(|x| COMP[*x as usize])
            .collect::<Vec<u8>>();
        Self {
            seq: revcomp,
            kmer_size: kmer_size as usize,
            start: 0,
        }
    }
}

impl Iterator for Reverse {
    type Item = Vec<u8>;

    /// Create a new kmer
    fn next(&mut self) -> Option<Self::Item> {
        if self.start + self.kmer_size > self.seq.len() {
            None
        } else {
            let kmer = self.seq[self.start..self.start + self.kmer_size].to_vec();
            self.start += 1;
            Some(kmer)
        }
    }
}

/// Struct that perform tokenization in canonical mode
struct Canonical<'a> {
    forward: Forward<'a>,
    reverse: Reverse,
}

impl<'a> Canonical<'a> {
    pub fn new(seq: &'a [u8], kmer_size: u64) -> Self {
        Self {
            forward: Forward::new(seq, kmer_size),
            reverse: Reverse::new(seq, kmer_size),
        }
    }
}

impl<'a> Iterator for Canonical<'a> {
    type Item = Vec<u8>;

    /// Create a new kmer
    fn next(&mut self) -> Option<Self::Item> {
        let forward = self.forward.next();
        let reverse = self.reverse.next();

        if forward < reverse {
            forward
        } else {
            reverse
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complement() {
        let mut lookup: [u8; 256] = [0; 256];

        for (i, item) in lookup.iter_mut().enumerate() {
            *item = i as u8;
        }

        lookup[b'A' as usize] = b'T';
        lookup[b'T' as usize] = b'A';
        lookup[b'C' as usize] = b'G';
        lookup[b'G' as usize] = b'C';
        lookup[b'a' as usize] = b'T';
        lookup[b't' as usize] = b'A';
        lookup[b'c' as usize] = b'G';
        lookup[b'g' as usize] = b'C';
        lookup[b'n' as usize] = b'N';

        assert_eq!(lookup, comp_lookup_table());
    }

    #[test]
    fn forward() {
        let tokenizer = Forward::new(b"AAGCTCAAGG", 5);

        assert_eq!(
            tokenizer.collect::<Vec<Vec<u8>>>(),
            vec![
                b"AAGCT".to_vec(),
                b"AGCTC".to_vec(),
                b"GCTCA".to_vec(),
                b"CTCAA".to_vec(),
                b"TCAAG".to_vec(),
                b"CAAGG".to_vec(),
            ]
        );
    }

    #[test]
    fn reverse() {
        let tokenizer = Reverse::new(b"AAGCTCAAGG", 5);

        assert_eq!(
            tokenizer.collect::<Vec<Vec<u8>>>(),
            vec![
                b"CCTTG".to_vec(),
                b"CTTGA".to_vec(),
                b"TTGAG".to_vec(),
                b"TGAGC".to_vec(),
                b"GAGCT".to_vec(),
                b"AGCTT".to_vec(),
            ]
        );
    }

    #[test]
    fn cannonical() {
        let tokenizer = Canonical::new(b"AAGCTCAAGG", 5);

        assert_eq!(
            tokenizer.collect::<Vec<Vec<u8>>>(),
            vec![
                b"AAGCT".to_vec(),
                b"AGCTC".to_vec(),
                b"GCTCA".to_vec(),
                b"CTCAA".to_vec(),
                b"GAGCT".to_vec(),
                b"AGCTT".to_vec(),
            ]
        );
    }

    #[test]
    fn tokenizer() -> error::Result<()> {
        let tokenizer = get_tokenizer(b"AAGCTCAAGG", 5, false, false)?;

        assert_eq!(
            tokenizer.collect::<Vec<Vec<u8>>>(),
            vec![
                b"AAGCT".to_vec(),
                b"AGCTC".to_vec(),
                b"GCTCA".to_vec(),
                b"CTCAA".to_vec(),
                b"GAGCT".to_vec(),
                b"AGCTT".to_vec(),
            ]
        );

        let tokenizer = get_tokenizer(b"AAGCTCAAGG", 5, true, false)?;

        assert_eq!(
            tokenizer.collect::<Vec<Vec<u8>>>(),
            vec![
                b"AAGCT".to_vec(),
                b"AGCTC".to_vec(),
                b"GCTCA".to_vec(),
                b"CTCAA".to_vec(),
                b"TCAAG".to_vec(),
                b"CAAGG".to_vec(),
            ]
        );

        let tokenizer = get_tokenizer(b"AAGCTCAAGG", 5, true, true)?;

        assert_eq!(
            tokenizer.collect::<Vec<Vec<u8>>>(),
            vec![
                b"CCTTG".to_vec(),
                b"CTTGA".to_vec(),
                b"TTGAG".to_vec(),
                b"TGAGC".to_vec(),
                b"GAGCT".to_vec(),
                b"AGCTT".to_vec(),
            ]
        );

        match get_tokenizer(b"AAGCTCAAGG", 5, false, true) {
            Err(ref e) if e.is::<error::Error>() => match e.downcast_ref::<error::Error>() {
                Some(error::Error::QueryReverseWithoutStranded) => (),
                _ => {
                    panic!("get_tokenizer should generate an error::Error::QueryReverseWithoutStranded")
                }
            },
            _ => {
                panic!("get_tokenizer should generate an error::Error::QueryReverseWithoutStranded")
            }
        }

        Ok(())
    }
}
