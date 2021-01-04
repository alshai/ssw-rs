use std::process;

struct SWOptions {
    a: i8,
    x: i8,
    o: i8,
    e: i8,
}

impl Default for SWOptions {
    fn default() -> Self {
        SWOptions {
            a: -1,
            x: 1,
            o: 1,
            e: 1,
        }
    }
}
// local alignment Smith Waterman algorithm
fn sw(opts: &SWOptions, t: &[u8], q: &[u8]) -> i8 {
    // closure for getting the index into 2D matrix given a row and column
    let nr = q.len() + 1;
    let nc = t.len() + 1;
    let get_idx = |r, c| r * nc + c;
    // initialize H, E and F matrices
    let mut hmat: Vec<i8> = vec![0; nc * nr];
    let mut emat: Vec<i8> = vec![0; nc * nr];
    let mut fmat: Vec<i8> = vec![0; nc * nr];
    for j in 1..nc { // initialize first row
        let idx = get_idx(0, j);
        hmat[idx] = 0;
        fmat[idx] = 0;
        emat[idx] = 0;
    }
    for i in 1..nr { // initialize first column
        let idx = get_idx(i, 0);
        hmat[idx] = 0;
        fmat[idx] = 0;
        emat[idx] = 0;
    }
    for i in 1..nr {
        for j in 1..nc {
            let idx = get_idx(i, j);
            let hscore = match q[i-1] == t[j-1] {
                true => hmat[get_idx(i-1, j-1)] + opts.a,
                false => hmat[get_idx(i-1, j-1)] + opts.x
            };
            emat[idx] = std::cmp::min(emat[get_idx(i-1,j)] + opts.e, hmat[get_idx(i-1,j)] + opts.o);
            fmat[idx] = std::cmp::min(fmat[get_idx(i,j-1)] + opts.e, hmat[get_idx(i,j-1)] + opts.o);
            let scores : [i8; 3] = [ fmat[idx], emat[idx], hscore ];
            hmat[idx] = *scores.iter().min().unwrap();
        }
    }
    return hmat[get_idx(nr-1, nc-1)];
}


pub fn sw_main() {
    let mut line = String::new(); 
    match std::io::stdin().read_line(&mut line) {
        Ok(_n) => {
            // split string into two
            let mut iter = line.split_whitespace();
            let str1 = iter.next().unwrap_or("");
            let str2 = iter.next().unwrap_or("");
            if str1 == "" || str2 == "" {
                process::exit(1);
            } 

            let sw_options: SWOptions = Default::default();
            let score = sw(&sw_options, str1.as_bytes(), str2.as_bytes());
            println!("{}", score);
        }
        Err(error) => {
            println!("error: {}", error);
        }
    }
}