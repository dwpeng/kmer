use kmer::*;

fn main() {
    let mut kc1 = KmerContainer::new(35);
    let seq = "ACTGATGCATGCTATCATCTACTATCATACTGTCTAGCTATCTATCCTTAGCTATATCA";
    kc1.eat(seq).unwrap();
    let mut ek = kc1.get_editable_kmer(0).unwrap();
    let kmer: Kmer = Kmer::from(&ek);
    println!("{}", kmer);
    println!("***************************");
    while ek.next().is_ok() {
        println!("{}", ek);
    }
}
