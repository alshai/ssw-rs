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
            a: -1,
            x: 1,
            o: 1,
            e: 1,
        }
    }
}
// local alignment Smith Waterman algorithm
#[allow(dead_code)]
pub fn sw(opts: &SWOptions, t: &[u8], q: &[u8]) -> i8 {
    // closure for getting the index into 2D matrix given a row and column
    let nr = q.len() + 1;
    let nc = t.len() + 1;
    let get_idx = |r, c| r * nc + c;
    // initialize H, E and F matrices
    let mut v_h: Vec<i8> = vec![0; nc * nr];
    let mut v_e: Vec<i8> = vec![0; nc * nr];
    let mut v_f: Vec<i8> = vec![0; nc * nr];
    for j in 1..nc { // initialize first row
        let idx = get_idx(0, j);
        v_h[idx] = 0;
        v_f[idx] = 0;
        v_e[idx] = 0;
    }
    for i in 1..nr { // initialize first column
        let idx = get_idx(i, 0);
        v_h[idx] = 0;
        v_f[idx] = 0;
        v_e[idx] = 0;
    }
    for i in 1..nr {
        for j in 1..nc {
            let idx = get_idx(i, j);
            let hscore = match q[i-1] == t[j-1] {
                true => v_h[get_idx(i-1, j-1)] + opts.a,
                false => v_h[get_idx(i-1, j-1)] + opts.x
            };
            v_e[idx] = std::cmp::min(v_e[get_idx(i-1,j)] + opts.e as i8, v_h[get_idx(i-1,j)] + opts.o as i8);
            v_f[idx] = std::cmp::min(v_f[get_idx(i,j-1)] + opts.e as i8, v_h[get_idx(i,j-1)] + opts.o as i8);
            let scores : [i8; 3] = [ v_f[idx], v_e[idx], hscore ];
            v_h[idx] = *scores.iter().min().unwrap();
        }
    }
    return v_h[get_idx(nr-1, nc-1)];
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
            for k in 0..15 { // loop over element in the register
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
    let x = std::mem::transmute::<__m128i, [u8; 16]>(_mm_cmpgt_epi8(a, b));
    for i in 0..x.len() {
        if x[i] > 0 {
            return true;
        }
    }
    return false;
}

#[allow(dead_code)]
// https://doi.org/10.1093/bioinformatics/btl582
pub unsafe fn ssw_8bit(opts: &SWOptions, t: &[u8], q: &[u8]) -> i8 {
    let _q = acgt_to_num(q);
    let _t = acgt_to_num(t);
    let v_profile = get_profile(&_q, opts.a, opts.x); // rows=chars, cols=seglen, elems=registers
    let seglen = ((_q.len()) + 15) / 16;
    println!("{} {}", _q.len(), seglen);
    let mut v_h_load =  vec![_mm_setzero_si128(); seglen];
    let mut v_h_store = vec![_mm_setzero_si128(); seglen];
    let mut v_max = _mm_setzero_si128();
    let mut v_e = vec![_mm_setzero_si128(); seglen];
    let v_go = _mm_set1_epi8(opts.o as i8);
    let v_ge = _mm_set1_epi8(opts.e as i8);
    for i in 0.._t.len() {
        let c = _t[i] as usize;
        let mut v_f = _mm_setzero_si128();
        let mut v_h = v_h_store[seglen-1];
        v_h = _mm_slli_si128(v_h, 1);
        // swap store and load
        let tmp = v_h_load;
        v_h_load = v_h_store;
        v_h_store = tmp;
        for j in 0..seglen {
            v_h = _mm_add_epi8(v_h, v_profile[c * seglen + j]);
            v_max = _mm_max_epi8(v_max, v_h);
            v_h = _mm_max_epi8(v_h, v_e[j]);
            v_h = _mm_max_epi8(v_h, v_f);
            v_h_store[j] = v_h;
            // calculate gap open/extend
            v_h = _mm_subs_epu8(v_h, v_go);
            v_e[j] = _mm_subs_epu8(v_e[j], v_ge);                
            v_e[j] = _mm_max_epi8(v_e[j], v_h);
            v_f = _mm_subs_epu8(v_f, v_ge);
            v_f = _mm_max_epi8(v_f, v_h);

            v_h = v_h_load[j];
        }
        // lazy f loop
        v_f = _mm_slli_si128(v_f, 1);
        let mut j = 0;
        while _mm_cmpgt_epi8_bool(v_f, _mm_subs_epu8(v_h_store[j], v_go)) {
            v_h_store[j] = _mm_max_epi8(v_h_store[j], v_f);
            v_f = _mm_subs_epu8(v_f, v_ge);
            j = j + 1;
            if j >= seglen {
                v_f = _mm_slli_si128(v_f, 1);
                j = 0;
            }
        }

    }
    return 0;
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