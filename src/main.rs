use kmer::*;

fn main() {
    let mut kc1 = KmerContainer::new(35);
    let seq = "ACTGATGCATGCTATCATCTACTATCATACTGTCTAGCTATCTATCCTTAGCTATATCA";
    kc1.eat(seq).unwrap();
    let mut ek = kc1.get_editable_kmer(0).unwrap();
    while ek.next().is_ok() {
        println!("{}", ek);
    }
}
