use std::{fmt::Display, vec};
use thiserror::Error;

use std::collections::VecDeque;

#[derive(Error, Debug)]
pub enum KmerError {
    #[error("Out of Index")]
    OutOfIndex,
    #[error("Too short sequence")]
    TooShortSequence,
}

#[derive(Debug)]
pub struct KmerContainer {
    pub __kmers: Vec<u8>,
    pub __block_size_per_kmer: usize,
    pub ksize: usize,
    pub size: usize,
    pub capacity: usize,
    // __phantom: std::marker::PhantomData<[u8; K]>,
}

static CAPACITY: usize = 1 << 5;
static BASE_MAP_NUMBER: [u8; 256] = {
    let mut map = [4; 256];
    map[b'A' as usize] = 0;
    map[b'C' as usize] = 1;
    map[b'G' as usize] = 2;
    map[b'T' as usize] = 3;
    map[b'a' as usize] = 0;
    map[b'c' as usize] = 1;
    map[b'g' as usize] = 2;
    map[b't' as usize] = 3;
    map
};
static NUMBER_MAP_BASE: [u8; 4] = [b'A', b'C', b'G', b'T'];

pub struct KmerBuffer<'a> {
    pub buff: VecDeque<u8>,
    pub sequence: &'a str,
    pub offset: usize,
    pub __ksize: usize,
    pub __kmer: Vec<u8>,
}

impl<'a> KmerBuffer<'a> {
    pub fn new(seq: &'a str, ksize: usize) -> Self {
        let mut buff = VecDeque::new();
        buff.push_back(0);
        for i in 0..ksize - 1 {
            buff.push_back(BASE_MAP_NUMBER[seq.as_bytes()[i] as usize]);
        }
        KmerBuffer {
            buff,
            sequence: seq,
            offset: ksize - 1,
            __ksize: ksize,
            __kmer: vec![0; ksize],
        }
    }

    pub fn reset(&mut self) {
        self.offset = self.__ksize - 1;
        self.buff.clear();
        self.buff.push_back(0);
        for i in 0..self.__ksize - 1 {
            self.buff
                .push_back(BASE_MAP_NUMBER[self.sequence.as_bytes()[i] as usize]);
        }
    }

    pub fn next(&mut self) -> Option<&[u8]> {
        if self.offset == self.sequence.len() {
            return None;
        }
        self.buff.pop_front();
        self.buff
            .push_back(BASE_MAP_NUMBER[self.sequence.as_bytes()[self.offset] as usize]);
        let mut kmer_ptr = self.__kmer.as_mut_ptr();
        for i in 0..self.__ksize {
            unsafe {
                *kmer_ptr = self.buff[i];
                kmer_ptr = kmer_ptr.add(1);
            };
        }
        self.offset += 1;
        Some(&self.__kmer)
    }
}

impl KmerContainer {
    pub fn new(ksize: usize) -> Self {
        let block_size = (ksize + 3) / 4;
        let kmers = vec![0; CAPACITY * block_size];
        KmerContainer {
            __kmers: kmers,
            __block_size_per_kmer: block_size,
            ksize: ksize,
            size: 0,
            capacity: CAPACITY,
        }
    }

    pub fn eat(&mut self, sequence: &str) -> Result<(), KmerError> {
        let seqlen = sequence.len();
        if seqlen < self.ksize {
            return Err(KmerError::TooShortSequence);
        }
        let nkmer = seqlen - self.ksize + 1;
        if self.size + nkmer < self.capacity {
            let new_capacity = self.capacity * 2;
            let new_size = new_capacity * self.__block_size_per_kmer;
            self.__kmers.resize(new_size, 0);
        }

        let mut kmer_ptr = self.__kmers.as_mut_ptr() as *mut u8;
        kmer_ptr = unsafe { kmer_ptr.add(self.size * self.__block_size_per_kmer) };
        let mut kmer_buff = KmerBuffer::new(sequence, self.ksize);
        for _ in 0..nkmer {
            let kmer = kmer_buff.next().unwrap();
            // pack kmer, every 4 bases in one u8
            // if less than 4 bases, fill with 0
            for i in 0..self.__block_size_per_kmer {
                unsafe {
                    *kmer_ptr = 0;
                    let mut pack = 0;
                    for j in 0..4 {
                        if i * 4 + j < self.ksize {
                            pack |= kmer[i * 4 + j] << (j * 2);
                        } else {
                            break;
                        }
                    }
                    *kmer_ptr = pack;
                    kmer_ptr = kmer_ptr.add(1);
                }
            }
            self.size += 1;
        }
        return Ok(());
    }

    pub fn get_editable_kmer(&mut self, index: usize) -> Result<EditableKmer, KmerError> {
        if index >= self.size {
            return Err(KmerError::OutOfIndex);
        }
        Ok(EditableKmer::from_kmer_container(self, index))
    }

    pub fn get_kmer(&self, index: usize) -> Result<Option<Kmer>, KmerError> {
        if index >= self.size {
            return Err(KmerError::OutOfIndex);
        }
        Ok(Kmer::from_kmer_container(self, index))
    }
}

#[derive(Debug)]
pub struct EditableKmer<'a> {
    pub __kc: &'a mut KmerContainer,
    pub __index: usize,
}

impl<'a> EditableKmer<'a> {
    pub fn from_kmer_container(kc: &'a mut KmerContainer, index: usize) -> Self {
        EditableKmer {
            __kc: kc,
            __index: index,
        }
    }

    #[inline]
    pub fn prev(&mut self) -> Result<(), KmerError> {
        if self.__index == 0 {
            return Err(KmerError::OutOfIndex);
        }
        self.__index -= 1;
        Ok(())
    }

    #[inline]
    pub fn next(&mut self) -> Result<(), KmerError> {
        if self.__index == self.__kc.size - 1 {
            return Err(KmerError::OutOfIndex);
        }
        self.__index += 1;
        Ok(())
    }

    #[inline]
    pub fn get_kmer_slice(&self) -> &[u8] {
        &self.__kc.__kmers[self.__index..self.__index + self.__kc.__block_size_per_kmer]
    }

    #[inline]
    pub fn get_kmer(&self) -> Kmer {
        Kmer::from_editable_kmer(self)
    }

    #[inline]
    pub fn ksize(&self) -> usize {
        self.__kc.ksize
    }
}

impl<'a> Display for EditableKmer<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let kmer = self.get_kmer();
        write!(f, "{}", kmer)
    }
}

#[derive(Debug)]
pub struct Kmer {
    pub bases: Vec<u8>,
    pub ksize: usize,
}

impl Display for Kmer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // every u8 have 4 bases
        for i in 0..self.ksize {
            write!(f, "{}", NUMBER_MAP_BASE[self.bases[i] as usize] as char)?;
        }
        Ok(())
    }
}

impl Kmer {
    pub fn from_editable_kmer(ek: &EditableKmer) -> Self {
        Kmer::from_kmer_container(ek.__kc, ek.__index).unwrap()
    }

    pub fn from_kmer_container(kc: &KmerContainer, index: usize) -> Option<Kmer> {
        if kc.size <= index {
            return None;
        }
        let mut bases = vec![0; kc.ksize];
        let nblock = (kc.ksize + 3) / 4;
        let offset = index * kc.__block_size_per_kmer;
        for i in 0..nblock {
            let mut packed = kc.__kmers[offset + i];
            for j in 0..4 {
                if i * 4 + j < kc.ksize {
                    bases[i * 4 + j] = packed & 0b11;
                    packed >>= 2;
                }
            }
        }
        Some(Kmer {
            bases,
            ksize: kc.ksize,
        })
    }
}

#[cfg(test)]
mod test {

    use super::*;

    #[test]
    fn test_new() {
        let kc = KmerContainer::new(3);
        assert_eq!(kc.size, 0);
        assert_eq!(kc.capacity, CAPACITY);
        assert_eq!(kc.__block_size_per_kmer, 1);

        let kc = KmerContainer::new(4);
        assert_eq!(kc.size, 0);
        assert_eq!(kc.capacity, CAPACITY);
        assert_eq!(kc.__block_size_per_kmer, 1);

        let kc = KmerContainer::new(5);
        assert_eq!(kc.size, 0);
        assert_eq!(kc.capacity, CAPACITY);
        assert_eq!(kc.__block_size_per_kmer, 2);
    }

    #[test]
    fn test_push() {
        let mut kc = KmerContainer::new(3);
        kc.eat("ACGT").unwrap();
        assert_eq!(kc.size, 2);
        assert!(kc.__kmers[0] == 0b100100);
        assert!(kc.__kmers[1] == 0b111001);
    }

    #[test]
    fn test_kmer() {
        let mut kc = KmerContainer::new(3);
        kc.eat("ACGT").unwrap();
        let kmer = kc.get_kmer(0).unwrap().unwrap();
        assert_eq!(format!("{}", kmer), "ACG");
        let kmer = kc.get_kmer(1).unwrap().unwrap();
        assert_eq!(format!("{}", kmer), "CGT");

        let mut kc = KmerContainer::new(4);
        kc.eat("ACGT").unwrap();
        let kmer = kc.get_kmer(0).unwrap().unwrap();
        assert_eq!(format!("{}", kmer), "ACGT");
        let kmer = kc.get_kmer(1);
        assert!(kmer.is_err());
        let mut kc = KmerContainer::new(11);
        let seq = "ACGTACGTACGtctactactagtctatct";
        kc.eat(seq).unwrap();
        assert!(kc.size == seq.len() - 11 + 1);
    }

    #[test]
    fn test_editable_kmer() {
        let mut kc = KmerContainer::new(11);
        let seq = "ACGTACGTACGtctactactagtctatct";
        kc.eat(seq).unwrap();

        let mut ek = kc.get_editable_kmer(0).unwrap();
        assert_eq!(format!("{}", ek), "ACGTACGTACG");

        assert!(ek.next().is_ok());
        let kmer = ek.get_kmer();
        assert_eq!(format!("{}", kmer), "CGTACGTACGT");
        assert_eq!(format!("{}", ek), "CGTACGTACGT");

        assert!(ek.prev().is_ok());
        assert_eq!(format!("{}", ek), "ACGTACGTACG");
    }
}
