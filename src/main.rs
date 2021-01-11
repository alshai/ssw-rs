mod sw;
mod tmp;

#[allow(dead_code)]
fn main() {
    // tmp::tmp();
    let opts = sw::SWOptions::default();
    let t = "AAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAA";
    let q = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
    println!("{}", t);
    println!("{}", q);
    let sw_score = sw::sw(&opts, t.as_bytes(), q.as_bytes());
    unsafe {
        let ssw_score = sw::ssw_8bit(&opts, t.as_bytes(), q.as_bytes());
        println!("sw: {}, ssw: {}", sw_score, ssw_score);
    }
}
