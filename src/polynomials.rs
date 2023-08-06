use blstrs::Scalar;
use ff::Field;
use std::ops::{AddAssign, Mul, MulAssign, SubAssign};
use more_asserts::debug_assert_le;
use crate::evaluation_domain::{BatchEvaluationDomain, EvaluationDomain};
use crate::{fft, is_power_of_two};

/// Returns the size of the evaluation domain needed to multiply these two polynomials via FFT.
pub fn get_evaluation_dom_size_for_multiplication(f: &Vec<Scalar>, g: &Vec<Scalar>) -> usize {
    let f_deg = f.len() - 1;
    let g_deg = g.len() - 1;

    // The degree $d$ of $f \cdot g$ will be $\deg{f} + \deg{g}$.
    let fg_deg = f_deg + g_deg;
    // But we need $d+1$ evaluations to interpolate a degree $d$ polynomial.
    let num_evals = fg_deg + 1;

    num_evals
}

/// Returns an evaluation domain for an FFT of size the number of coefficients in the polynomial $f(X) \cdot g(X)$.
pub fn get_evaluation_dom_for_multiplication(f: &Vec<Scalar>, g: &Vec<Scalar>) -> EvaluationDomain {
    let num_evals = get_evaluation_dom_size_for_multiplication(f, g);
    EvaluationDomain::new( num_evals).unwrap()
}

/// Returns $f(x)$ for a polynomial $f$ and a point $x$.
pub fn poly_eval(f: &Vec<Scalar>, x: &Scalar) -> Scalar {
    assert!(!f.is_empty());

    //let deg = f.len() - 1;
    let mut eval = Scalar::zero();  // f(x)
    let mut x_i = Scalar::one();    // x^i, i = {0, 1, ..., deg(f)}
    for c_i in f {
        eval += c_i * x_i;

        x_i *= x;
    }

    eval
}

/// Lets $f(X) = f(X) + g(X)$. Stores the result in `f` so assumes $\deg(f) \ge \deg(g)$.
pub fn poly_add_assign(f: &mut Vec<Scalar>, g: &[Scalar]) {
    if g.len() > f.len() {
        f.resize(g.len(), Scalar::zero());
    }

    for i in 0..g.len() {
        f[i].add_assign(g[i]);
    }
}

/// Lets $f(X) = f(X) - g(X)$. Stores the result in `f` so assumes $\deg(f) \ge \deg(g)$.
pub fn poly_sub_assign(f: &mut Vec<Scalar>, g: &[Scalar]) {
    if g.len() > f.len() {
        f.resize(g.len(), Scalar::zero());
    }
    for i in 0..g.len() {
        f[i].sub_assign(g[i]);
    }
}

/// Computes the product of $f$ and $g$, letting $f = f \cdot g$ and $g = FFT(g)$.
/// Let $d = \deg(f) + \deg(g)$. Takes $O(d\log{d})$ time via three FFT.
///
/// Note: If you already have an `EvaluationDomain` for $n = \deg(f) + \deg(g) + 1$, use the faster
/// `poly_mul_assign_fft_with_dom` which avoids recomputing the roots of unity and the other inverses.
pub fn poly_mul_assign_fft(f: &mut Vec<Scalar>, g: &mut Vec<Scalar>) {
    debug_assert!(!f.is_empty());
    debug_assert!(!g.is_empty());

    poly_mul_assign_fft_with_dom(f, g, &get_evaluation_dom_for_multiplication(f, g))
}

/// Like `poly_mul_assign` but returns the result, instead of modifying the arguments.
pub fn poly_mul_fft(f: &Vec<Scalar>, g: &Vec<Scalar>) -> Vec<Scalar> {
    debug_assert!(!f.is_empty());
    debug_assert!(!g.is_empty());

    let mut f_copy = f.clone();
    let mut g_copy = g.clone();

    poly_mul_assign_fft_with_dom(&mut f_copy, &mut g_copy, &get_evaluation_dom_for_multiplication(f, g));

    f_copy
}

/// If the caller already has an `EvaluationDomain` for $n = \deg(f) + \deg(g) + 1$, this function
/// will avoid some redundant field operations and be slighlty faster than `poly_mul_assign_fft`.
pub fn poly_mul_assign_fft_with_dom(f: &mut Vec<Scalar>, g: &mut Vec<Scalar>, dom: &EvaluationDomain) {
    debug_assert!(!f.is_empty());
    debug_assert!(!g.is_empty());
    debug_assert_eq!((f.len() - 1) + (g.len() - 1) + 1, dom.n);

    fft::fft_assign(f, &dom);
    fft::fft_assign(g, &dom);
    for i in 0..dom.N {
        f[i].mul_assign(g[i]);
    }

    fft::ifft_assign(f, &dom);
    f.truncate(dom.n);
}

/// Like `poly_mul_assign_fft` but slower in time $\deg(f) \cdot \deg(g)$ and returns the product in `out`, leaving `f` and `g` untouched.
/// TODO(Perf): Not sure if we can do this in-place over `f` or `g` without a separate `out`.
pub fn poly_mul_assign_slow(f: &Vec<Scalar>, g: &Vec<Scalar>, out: &mut Vec<Scalar>)  {
    assert!(!f.is_empty());
    assert!(!g.is_empty());

    let f_len = f.len();
    let g_len = g.len();

    // Set `out` to all zeros.
    out.truncate(0);    // Rust docs say "Note that this method has no effect on the allocated capacity of the vector."
    out.resize(f_len + g_len - 1, Scalar::zero());
    for (i1, a) in f.iter().enumerate() {
        for (i2, b) in g.iter().enumerate() {
            let mut prod = *a;
            prod.mul_assign(b);
            out[i1 + i2].add_assign(&prod);
        }
    }
}

/// Like `poly_mul_fft` but slower: runs in time $\deg(f) \cdot \deg(g)$.
/// Returns the product, leaving `f` and `g` untouched.
///
/// Performance notes: Beats FFT for $t \ge 32$.
pub fn poly_mul_slow(f: &Vec<Scalar>, g: &Vec<Scalar>) -> Vec<Scalar> {
    let mut out = Vec::with_capacity(get_evaluation_dom_size_for_multiplication(f, g));
    poly_mul_assign_slow(f, g, &mut out);
    out
}

/// Like `poly_mul_assign_fft` but runs in time $O(n)^{1.58})$.
/// Assumes that $\deg{f} - 1 = \deg{g} - 1 = 2^k$ for some $k$.
/// Returns the product in `out`, leaving `f` and `g` untouched.
///
/// Performance notes: Starts beating `poly_mul_slow` after $t \ge 512$. Never beats FFT, so useless
/// because FFT beats naive after $t > 64$.
pub fn poly_mul_assign_less_slow(f: &Vec<Scalar>, g: &Vec<Scalar>, out: &mut Vec<Scalar>) {
    let mut result = poly_mul_less_slow(f.as_slice(), g.as_slice());

    out.truncate(0);
    out.append(&mut result);
}

pub fn poly_mul_less_slow(f: &[Scalar], g: &[Scalar]) -> Vec<Scalar> {
    let n = f.len();
    assert!(is_power_of_two(n));

    // TODO(Perf): Let the base case be the naive n^2 multiplication for any f_len < threshold?
    // handle base case
    if n == 1 {
        return vec![f[0].mul(g[0])];
    }

    let g_len = g.len();
    let half_n = n / 2;
    assert_eq!(n, g_len);

    // split f into halves
    let f_0 = &f[0..half_n];
    let f_1 = &f[half_n..n];
    debug_assert_eq!(f_0.len() + f_1.len(), n);

    // split g into halves
    let g_0 = &g[0..half_n];
    let g_1 = &g[half_n..n];
    debug_assert_eq!(g_0.len() + g_1.len(), n);

    // let f_01 = f_0 + f_1
    let mut f_01 = Vec::with_capacity(half_n);
    f_0.clone_into(&mut f_01);
    for i in 0..half_n {
        f_01[i].add_assign(f_1[i]);
    }

    // let g_01 = g_0 + g_1
    let mut g_01 = Vec::with_capacity(half_n);
    g_0.clone_into(&mut g_01);
    for i in 0..half_n {
        g_01[i].add_assign(g_1[i]);
    }

    // $y = (f_0 + f_1) \cdot (g_0 + g_1)$
    let mut y = poly_mul_less_slow(&f_01, &g_01);
    // $z = f_0 \cdot g_0$
    let mut z = poly_mul_less_slow(f_1, g_1);
    // $u = f_1 \cdot g_0$
    let mut u = poly_mul_less_slow(f_0, g_0);

    // first, compute (y - u - z)
    poly_sub_assign(&mut y, &u);
    poly_sub_assign(&mut y, &z);

    // second, compute (y - u - z) * X^{n/2}
    poly_xnmul_assign(&mut y, half_n);

    // third, compute z * X^n
    poly_xnmul_assign(&mut z, n);

    // fourth, add everything up with $u(X)$
    poly_add_assign(&mut u, &y);
    poly_add_assign(&mut u, &z);

    u.into()
}
/// Sets $f(X) = f(X) \cdot X^n$, by simply shifting the coefficients.
/// As always we assume $\deg{f}$ is `f.len() - 1`.
pub fn poly_xnmul_assign(f: &mut Vec<Scalar>, n: usize) {
    if n == 0 {
        return
    }

    let old_len = f.len();

    // extend with zero coefficients for X^n, X^{n-1}, \dots, X
    f.resize(old_len + n, Scalar::zero());

    // Shift coefficients by `n` positions
    for i in (0..old_len).rev() {
        f[i + n] = f[i]
    }

    // Set the last n coefficients $f_{n-1}, \cdots, f_0$ to 0.
    for i in 0..n {
        f[i] = Scalar::zero();
    }
}

/// Like `poly_mul_by_xn_assign` but returns the result.
pub fn poly_xnmul(f: &Vec<Scalar>, n: usize) -> Vec<Scalar> {
    if n == 0 {
        return f.clone();
    }

    let len = n + f.len();

    // let result = f
    let mut result = Vec::with_capacity(len);
    result.resize(f.len(), Scalar::zero());
    result.copy_from_slice(f);

    poly_xnmul_assign(&mut result, n);

    result
}

/// TODO(Perf): Multi-thread the tree-based computation.
/// TODO(Perf): Avoid allocations; do it in place. Trickier due to the 2x larger FFTs and the slow multiplication algorithm.
/// Given a set $S$ of scalars, returns the *accumulator* polynomial $Z(X) = \prod_{a \in S} (X - a)$.
#[allow(non_snake_case)]
pub fn accumulator_poly_slow(S: &[Scalar]) -> Vec<Scalar> {
    let set_size = S.len();
    // println!("len: {len}");

    if set_size == 0 {
        return vec![];
    } else if set_size == 1 {
        // println!("Returning (X - {})", S[0]);
        return vec![-S[0], Scalar::one()]
    }

    let m = set_size / 2;
    // println!("m: {m}");

    let left =  &S[0..m];
    // println!("left: {}", left.len());
    let right = &S[m..set_size];
    // println!("right: {}", right.len());

    let mut left_poly = accumulator_poly_slow(left);
    let mut right_poly = accumulator_poly_slow(right);
    if left_poly.is_empty() {
        return right_poly
    }
    if right_poly.is_empty() {
        return left_poly
    }

    poly_mul_assign_fft(&mut left_poly, &mut right_poly);

    left_poly
}

/// Like `accumulator_poly_slow` but:
///  - Avoids recomputing too many different roots of unity (EvaluationDomain::new takes 3.5 microsecs
///    and `accumulator_poly_slow` makes around 2000 calls to it). Saves 10 milliseconds out of total of 70.
///  - Avoids FFTs when multiplying polynomials with $\deg{f} + \deg{g} - 1 \le fft_thresh$.
///    Saves 35 milliseconds.
#[allow(non_snake_case)]
pub fn accumulator_poly(S: &[Scalar], batch_dom: &BatchEvaluationDomain, fft_thresh: usize) -> Vec<Scalar> {
    let set_size = S.len();

    if set_size == 0 {
        return vec![];
    } else if set_size == 1 {
        return vec![-S[0], Scalar::one()]
    } else if set_size == 2 {
        return vec![S[0]*S[1], -(S[0] + S[1]), Scalar::one()];
    } else if set_size == 3 {
        // https://math.stackexchange.com/questions/88917/relation-betwen-coefficients-and-roots-of-a-polynomial
        let s1_add_s2 = S[1] + S[2];
        let s1_mul_s2 = S[1] * S[2];

        let c_0 = -(S[0] * s1_mul_s2);  // -(S[0] S[1] S[2])
        let c_1 = S[0] * s1_add_s2 + s1_mul_s2; // S[0] S[1] + S[0] S[2] + S[1] S[2]
        let c_2 = -(S[0] + S[1] + S[2]);

        return vec![c_0, c_1, c_2, Scalar::one()];
    }

    let m = set_size / 2;

    let left =  &S[0..m];
    let right = &S[m..set_size];

    let mut left_poly = accumulator_poly(left, batch_dom, fft_thresh);
    let mut right_poly = accumulator_poly(right, batch_dom, fft_thresh);
    if left_poly.is_empty() {
        return right_poly
    }
    if right_poly.is_empty() {
        return left_poly
    }

    let dom_size = get_evaluation_dom_size_for_multiplication(&left_poly, &right_poly);
    if dom_size < fft_thresh {
        poly_mul_slow(&left_poly, &right_poly)
    } else {
        poly_mul_assign_fft_with_dom(&mut left_poly, &mut right_poly, &batch_dom.get_subdomain(dom_size));
        left_poly
    }
}

/// More wisely schedules the sizes of the polynomials that are multiplied via FFT. Specifically,
/// always ensures we are multiplying degree (2^k - 1) by degree 2^k, to get a degree (2^{k+1} - 1)
/// polynomial. This ensures optimal usage of the 3 FFTs: the first two FFTs on the degree 2^k - 1
/// and degree 2^k have 50% usage, while the last inverse FFT that outputs 2^{k+1} coefficients has
/// 100% usage.
///
/// Expected this to be much faster than `accumulator_poly` but seems to beat by at most 1 millisecond
/// for the following parameters:
///
///   128 FFT thresh, 256 naive -thresh > 15.4 ms
///   256 FFT thresh, 512 naive thresh -> 14.8 ms
///   256 FFT thresh, 128 naive thresh -> 14.1 ms
#[allow(non_snake_case)]
pub fn accumulator_poly_scheduled(S: &[Scalar], batch_dom: &BatchEvaluationDomain, naive_thresh: usize, fft_thresh: usize) -> Vec<Scalar> {
    let mut n = S.len() + 1;

    if S.len() < naive_thresh {
        println!("Directly returning accumulator for set size {}", S.len());
        return accumulator_poly(S, batch_dom, fft_thresh);
    }

    let mut batch_size = 1;
    while n != 1 {
        n /= 2;
        batch_size *= 2;
    }
    batch_size -= 1;    // 2^k - 1 <= |S|

    debug_assert_le!(batch_size, S.len());
    assert!(is_power_of_two(batch_size + 1));

    let left = accumulator_poly_scheduled_inner(&S[0..batch_size], batch_dom, naive_thresh, fft_thresh);
    if batch_size == S.len() {
        left
    } else {
        let right = accumulator_poly_scheduled(&S[batch_size..], batch_dom, naive_thresh, fft_thresh);

        poly_mul_fft(&left, &right)
    }
}

/// Calculates the derivative of $f(X)$.
pub fn poly_differentiate(f: &mut Vec<Scalar>) {
    let f_deg = f.len() - 1;

    for i in 0..f_deg {
        f[i] = f[i+1].mul(Scalar::from((i + 1) as u64));
    }

    f.truncate(f_deg);
}

#[allow(non_snake_case)]
fn accumulator_poly_scheduled_inner(S: &[Scalar], batch_dom: &BatchEvaluationDomain, naive_thresh: usize, fft_thresh: usize) -> Vec<Scalar> {
    let len = S.len();
    debug_assert!(is_power_of_two(len + 1));

    if len < naive_thresh {
        println!("Directly returning accumulator for set size {len}");
        return accumulator_poly(S, batch_dom, fft_thresh);
    }

    let batch_size = (len + 1) / 2 - 1;
    debug_assert_eq!(batch_size * 2 + 1, len);

    let mut b1 = accumulator_poly_scheduled_inner(&S[0..batch_size], batch_dom, naive_thresh, fft_thresh);
    debug_assert_eq!(b1.len(), batch_size + 1);
    let b2 = accumulator_poly_scheduled_inner(&S[batch_size..2* batch_size], batch_dom, naive_thresh, fft_thresh);
    debug_assert_eq!(b2.len(), batch_size + 1);
    let deg1 = accumulator_poly(&S[2* batch_size..], batch_dom, fft_thresh);
    debug_assert_eq!(deg1.len(), 2);

    println!("Multiplying deg-{} by deg-{}", b1.len() - 1, b2.len() - 1 + deg1.len() -1);

    let mut b2 = poly_mul_slow(&deg1, &b2);
    poly_mul_assign_fft_with_dom(&mut b1, &mut b2, &batch_dom.get_subdomain(len + 1));

    b1
}