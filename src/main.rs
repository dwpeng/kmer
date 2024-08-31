use kmer::*;



fn main() {
    let kc1 = KmerContainer::<3>::new();
    println!("{:?}", kc1);
    let kmer1 = Kmer::<3>::from_kmer_container(&kc1, 0);
    println!("Hello, world! {:?}", kmer1)
}
