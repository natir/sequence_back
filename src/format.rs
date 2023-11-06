//! Detect format of file

/* std use */

/* crate use */

/* project use */
use crate::error;

/// List available reads format
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub enum ReadsFormat {
    /// Fasta format
    #[default]
    Fasta,
    /// Fastq format
    Fastq,
}

impl std::fmt::Display for ReadsFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl ReadsFormat {
    /// Detect format of file by read the first byte
    pub fn detect<R>(input: &mut R) -> error::Result<Self>
    where
        R: std::io::Read,
    {
        let mut first_byte: [u8; 1] = [0];

        input.read_exact(&mut first_byte)?;

        match first_byte {
            [b'>'] => Ok(Self::Fasta),
            [b'@'] => Ok(Self::Fastq),
            _ => Err(error::Error::NotFastaOrFastq.into()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detect_fasta() -> error::Result<()> {
        let mut input = std::io::Cursor::new([b'>']);

        assert_eq!(ReadsFormat::detect(&mut input)?, ReadsFormat::Fasta);

        Ok(())
    }

    #[test]
    fn detect_fastq() -> error::Result<()> {
        let mut input = std::io::Cursor::new([b'@']);

        assert_eq!(ReadsFormat::detect(&mut input)?, ReadsFormat::Fastq);

        Ok(())
    }

    #[test]
    fn read_error() -> error::Result<()> {
        let mut input = std::io::Cursor::new([]);

        match ReadsFormat::detect(&mut input) {
            Ok(_) => panic!("ReadsFormat::detect should generate an std::io::Error"),
            Err(e) if e.is::<std::io::Error>() => (),
            Err(_) => panic!("ReadsFormat::detect should generate an std::io::Error"),
        }

        Ok(())
    }

    #[test]
    fn format_error() -> error::Result<()> {
        let mut input = std::io::Cursor::new([b'A']);

        match ReadsFormat::detect(&mut input) {
            Err(ref e) if e.is::<error::Error>() => match e.downcast_ref::<error::Error>() {
                Some(error::Error::NotFastaOrFastq) => (),
                _ => {
                    panic!("ReadsFormat::detect should generate an error::Error::NotAFastaOrFastq")
                }
            },
            _ => {
                panic!("ReadsFormat::detect should generate an error::Error::NotAFastaOrFastq")
            }
        }

        Ok(())
    }

    #[test]
    fn display() {
        assert_eq!(format!("{}", ReadsFormat::Fasta), "Fasta");
        assert_eq!(format!("{}", ReadsFormat::Fastq), "Fastq");
    }
}
