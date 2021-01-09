mod sw;
mod tmp;

#[allow(dead_code)]
fn main() {
    // tmp::tmp();
    let opts = sw::SWOptions::default();
    unsafe {
        let t = "AAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAA";
        let q = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
        println!("{:?}", sw::acgt_to_num(t.as_bytes()));
        println!("{:?}", sw::acgt_to_num(q.as_bytes()));
        sw::ssw_8bit(&opts, t.as_bytes(), q.as_bytes());
    }
}
