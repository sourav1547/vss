use std::ops::Mul;
use blstrs::Scalar;
use ff::Field;
use rand::thread_rng;
use vss::{random_scalar, random_scalars};
use vss::evaluation_domain::{BatchEvaluationDomain, EvaluationDomain};
use vss::fft::{fft_assign, ifft_assign};
use vss::polynomials::{poly_eval, poly_mul_fft, poly_mul_less_slow, poly_mul_slow, poly_xnmul};

#[test]
#[allow(non_snake_case)]
fn test_batch_evaluation_domain() {
    for N in 1..=16 {
        let batch_dom = BatchEvaluationDomain::new(N);

        for k in 1..=N {
            let dom1 = EvaluationDomain::new(k).unwrap();
            let dom2 = batch_dom.get_subdomain(k);
            assert_eq!(dom1, dom2);
        }
    }
}

#[test]
#[allow(non_snake_case)]
fn test_fft_assign() {
    let mut rng = thread_rng();
    for N in [7, 8, 16] {
        let dom = EvaluationDomain::new(N).unwrap();

        let mut f = random_scalars(N, &mut rng);
        //println!("f(X): {:#?}", f);
        let f_orig = f.clone();

        // Computes $f(\omega^i)$ for all $i \in N$, where $\omega$ is an $N$th root of unity and
        // $N$ is the smallest power of two larger than $n$ (or $n$ itself if $n=2^k$).
        fft_assign(&mut f, &dom);
        //println!("FFT(f(X)): {:#?}", f);

        // Correctness test #1: Test against inverse FFT.
        let mut f_inv = f.clone();
        ifft_assign(&mut f_inv, &dom);
        //println!("FFT^{{-1}}(FFT(f(X))): {:#?}", f_inv);

        if f_inv.len() > f_orig.len() {
            let mut i = f_orig.len();
            while i < f_inv.len() {
                assert_eq!(f_inv[i], Scalar::zero());
                i += 1;
            }
            f_inv.truncate(f_orig.len());
        }
        assert_eq!(f_inv, f_orig);

        // Correctness test #2: Test by re-evaluating naively at the roots of unity
        let mut omega = Scalar::one();
        for i in 0..N {
            let y_i = poly_eval(&f_orig, &omega);

            //println!("y[{i}]: {}", y_i);
            assert_eq!(y_i, f[i]);

            omega *= dom.get_primitive_root_of_unity();
        }
    }
}

#[test]
fn test_poly_multiply() {
    let mut rng = thread_rng();
    for num_coeffs_f in [1, 2, 3, 4, 5, 6, 7, 8] {
        for num_coeffs_g in [1, 2, 3, 4, 5, 6, 7, 8] {
            let f = random_scalars(num_coeffs_f, &mut rng);
            let g = random_scalars(num_coeffs_g, &mut rng);

            // FFT-based multiplication
            let fft_fg = poly_mul_fft(&f, &g);

            // Naive multiplication
            let naive_fg = poly_mul_slow(&f, &g);

            // We test correctness of $h(X) = f(X) \cdot g(X)$ by picking a random point $r$ and
            // comparing $h(r)$ with $f(r) \cdot g(r)$.
            let r = random_scalar(&mut rng);

            let fg_rand = poly_eval(&f, &r).mul(poly_eval(&g, &r));
            let fft_fg_rand = poly_eval(&fft_fg, &r);
            assert_eq!(fft_fg_rand, fg_rand);

            // We also test correctness of the naive multiplication algorithm
            let naive_fg_rand = poly_eval(&naive_fg, &r);
            assert_eq!(naive_fg_rand, fg_rand);

            // Lastly, of course the naive result should be the same as the FFT result (since they are both correct)
            assert_eq!(naive_fg, fft_fg);
        }
    }
}

#[test]
fn test_poly_multiply_divide_and_conquer() {
    let mut rng = thread_rng();
    for log_n in [1, 2, 3, 4, 5, 6, 7, 8] {
        let n = 1 << log_n;
        let f = random_scalars(n, &mut rng);
        let g = random_scalars(n, &mut rng);

        let fg = poly_mul_less_slow(&f, &g);

        // FFT-based multiplication
        let fft_fg = poly_mul_fft(&f, &g);
        assert_eq!(fg, fft_fg);

        // Schwartz-Zippel test
        let r = random_scalar(&mut rng);
        let fg_rand = poly_eval(&f, &r).mul(poly_eval(&g, &r));
        let our_fg_rand = poly_eval(&fg, &r);
        assert_eq!(our_fg_rand, fg_rand);
    }
}

#[test]
#[allow(non_snake_case)]
fn test_poly_shift() {
    let mut rng = thread_rng();
    for num_coeffs_f in [1, 2, 3, 4, 5, 6, 7, 8] {
        for n in 0..16 {
            // compute the coefficients of X^n
            let mut Xn = Vec::with_capacity(n + 1);
            Xn.resize(n + 1, Scalar::zero());
            Xn[n] = Scalar::one();

            // pick a random f
            let f = random_scalars(num_coeffs_f, &mut rng);

            // f(X) * X^n via shift
            let shifted1 = poly_xnmul(&f, n);
            // f(X) * X^n via multiplication
            let shifted2 = poly_mul_fft(&f, &Xn);

            assert_eq!(shifted1, shifted2);
        }
    }
}