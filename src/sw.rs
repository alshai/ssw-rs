use std::process;


#[cfg(target_arch = "x86")]
use std::arch::x86::*;
#[allow(unused_imports)]
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

const E: u8 = 4;


pub struct SWOptions {
    a: i8,
    x: i8,
    o: u8,
    e: u8,
}

impl Default for SWOptions {
    fn default() -> Self {
        SWOptions {
            a: 1, // match (pos or neg)
            x: -1, // mismatch (pos or neg)
            o: 1, // gap open (must be pos)
            e: 1, // gap extend (must be pos)
        }
    }
}
// (unoptimized) local alignment Smith Waterman algorithm
pub fn sw(opts: &SWOptions, t: &[u8], q: &[u8]) -> i8 {
    // closure for getting the index into 2D matrix given a row and column
    let nr = q.len() + 1;
    let nc = t.len() + 1;
    let get_idx = |r, c| r * nc + c;
    // initialize H, E and F matrices
    let mut m_h: Vec<i8> = vec![0; nc * nr];
    let mut m_e: Vec<i8> = vec![0; nc * nr];
    let mut m_f: Vec<i8> = vec![0; nc * nr];
    let mut max_score = 0;
    for j in 1..nc { // initialize first row
        let idx = get_idx(0, j);
        m_h[idx] = 0;
        m_f[idx] = 0;
        m_e[idx] = 0;
    }
    for i in 1..nr { // initialize first column
        let idx = get_idx(i, 0);
        m_h[idx] = 0;
        m_f[idx] = 0;
        m_e[idx] = 0;
    }
    for i in 1..nr {
        for j in 1..nc {
            let idx = get_idx(i, j);
            let hscore = match q[i-1] == t[j-1] {
                true => m_h[get_idx(i-1, j-1)] + opts.a,
                false => m_h[get_idx(i-1, j-1)] + opts.x
            };
            m_e[idx] = std::cmp::max(m_e[get_idx(i-1,j)] - opts.e as i8, m_h[get_idx(i-1,j)] - opts.o as i8);
            m_f[idx] = std::cmp::max(m_f[get_idx(i,j-1)] - opts.e as i8, m_h[get_idx(i,j-1)] - opts.o as i8);
            let scores : [i8; 4] = [0, m_f[idx], m_e[idx], hscore ];
            m_h[idx] = *scores.iter().max().unwrap();
            max_score = std::cmp::max(m_h[idx], max_score);
        }
    }
    // for i in 0..nr {
    //     for j in 0..nc {
    //         print!("{} ", m_h[get_idx(i,j)]);
    //     }
    //     print!("\n");
    // }
    // TODO: do traceback :)
    return max_score;
}


#[allow(dead_code)]
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

// TODO: enable support for using custom alphabets and a substution matrix
// q: query sequence as numbers (acgt_to_num)
// a: match score
// x: mismatch score
#[allow(dead_code)]
pub unsafe fn get_profile(q: &[u8], a: i8, x: i8) -> std::vec::Vec<__m128i> {
    for i in 0..q.len() {
        if q[i] >= E { 
            panic!("get_profile error: expected query element to be less than E. Got _ instead. Did you use acgt_to_num() on query beforehand?");
        }
    }
    let seglen : usize = (q.len() + 15) / 16;
    // A: [16][16][16][16] ...
    // C: [16][16][16][16] ...
    // G: [16][16][16][16] ...
    // T: [16][16][16][16] ..seglen.. [15]
    // rows=ncharacter, cols = seglen
    let mut profile = vec![_mm_setzero_si128(); (E as usize) * seglen];
    for c in 0u8..E {
        for i in 0..seglen { // loop over columns (each col is register)
            let mut j = i; // for accessing query chars
            let regidx = seglen * (c as usize) + i;
            let mut reg = vec![0i8; 16];
            for k in 0..16 { // loop over element in the register
                if j >= q.len() {
                    reg[k] = 0;
                }
                else {
                    reg[k] = match q[j] == c {
                        true => a,
                        false => x, 
                    }
                }
                j = j + seglen; // ensure that we're following striped pattern along query (0, seglen, 2*seglen, etc)
            }
            // turn the vec into a register and add to profile
            profile[regidx] = _mm_load_si128(&reg[0] as *const _ as *const __m128i);
        }
    }
    return profile;
}

unsafe fn _mm_cmpgt_epi8_bool(a: __m128i, b: __m128i) -> bool {
    let x = std::mem::transmute::<__m128i, [i8; 16]>(_mm_cmpgt_epi8(a, b));
    for i in 0..x.len() {
        if x[i] > 0 {
            return true;
        }
    }
    return false;
}

#[allow(dead_code)]
unsafe fn print_8bitreg(r: __m128i) {
    let v = std::mem::transmute_copy::<__m128i, [i8; 16]>(&r);
    println!("{:?}", v);
}

// https://doi.org/10.1093/bioinformatics/btl582
pub unsafe fn ssw_8bit(opts: &SWOptions, t: &[u8], q: &[u8]) -> i8 {
    let _q = acgt_to_num(q);
    let _t = acgt_to_num(t);
    let v_profile = get_profile(&_q, opts.a, opts.x); // rows=chars, cols=seglen, elems=registers
    let seglen = ((_q.len()) + 15) / 16;
    let mut v_h_load =  vec![_mm_setzero_si128(); seglen];
    let mut v_h_store = vec![_mm_setzero_si128(); seglen];
    let mut v_max = _mm_set1_epi8(-128);
    let mut v_e = vec![_mm_setzero_si128(); seglen];
    let v_go = _mm_set1_epi8(opts.o as i8);
    let v_ge = _mm_set1_epi8(opts.e as i8);
    for i in 0.._t.len() {
        let c = _t[i] as usize;
        let mut v_f = _mm_setzero_si128();
        // align previous iter's bottom register to current iter's top register by shifting it down
        // ie. set v_h to equiv of v_h_load[-1] << 1
        let mut v_h = v_h_store[seglen-1]; 
        v_h = _mm_slli_si128(v_h, 1); // adj s.t. i==1 takes from 0
        // load and adjust upper-left diagonal from last vec of previous iter
        std::mem::swap(&mut v_h_load, &mut v_h_store);
        for j in 0..seglen { // populate v_h_store
            // update v_h for current current j
            v_h = _mm_adds_epi8(v_h, v_profile[c * seglen + j]); // calc match/mismatch
            // add gap penalties that we calculated last iteration
            v_h = _mm_max_epi8(v_h, v_e[j]); 
            v_h = _mm_max_epi8(v_h, v_f);
            v_h = _mm_max_epi8(v_h, _mm_setzero_si128()); // this is not in pseudocode bc we don't bias scores (we assume sse4)
            v_max = _mm_max_epi8(v_max, v_h); // keep track of overall max
            v_h_store[j] = v_h; // store v_h for next iteration
            // updated v_f and v_e (gap open/extend) for next iteration
            v_h = _mm_subs_epi8(v_h, v_go); // calulate open penalty
            v_e[j] = _mm_subs_epi8(v_e[j], v_ge); // insertion extend penalty
            v_e[j] = _mm_max_epi8(v_e[j], v_h); // determine insertion open v extend
            v_f = _mm_subs_epi8(v_f, v_ge); // deletion extend penalty
            v_f = _mm_max_epi8(v_f, v_h); // determine delete open v extend
            // get ready for next iteration (1 cell down)
            v_h = v_h_load[j]; // next iter will use v_h_load[j-1]
        }
        // lazy f loop - update v_h_store if v_f affects anything
        v_f = _mm_slli_si128(v_f, 1);
        let mut j = 0;
        while _mm_cmpgt_epi8_bool(v_f, _mm_subs_epi8(v_h_store[j], v_go)) {
            v_h_store[j] = _mm_max_epi8(v_h_store[j], v_f);
            v_f = _mm_subs_epi8(v_f, v_ge);
            j = j + 1;
            if j >= seglen {
                v_f = _mm_slli_si128(v_f, 1);
                j = 0;
            }
        }
    }
    // get overall max from v_max
    let v_max_vec = std::mem::transmute::<__m128i, [i8; 16]>(v_max);
    match v_max_vec.iter().max() {
        Some(max) => return *max,
        None => panic!("failed to to ssw_8bit"),
    }
}

#[allow(dead_code)]
pub fn acgt_to_num(s: &[u8]) -> std::vec::Vec<u8> {
    let mut v = vec![0; s.len()];
    for (i, c) in s.iter().enumerate() {
        v[i] = match c.to_ascii_uppercase() {
            65 => 0,
            67 => 1, 
            71 => 2, 
            84 => 3,
            _ => panic!("Non-ACGT value found"),
        }
    }
    return v;
}