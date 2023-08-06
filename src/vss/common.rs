use std::ops::Mul;

use blstrs::{Scalar, G1Projective};
use group::Group;
use rand::{distributions::Uniform, prelude::Distribution, thread_rng};

use crate::{evaluation_domain::BatchEvaluationDomain, lagrange::all_lagrange_denominators, random_scalars, fft::fft_assign, pvss::SharingConfiguration};


/// Return a random scalar within a small range [0,n) 
pub fn random_scalar_range<R>(mut rng: &mut R, u: u64) -> Scalar 
    where R: rand_core::RngCore + rand::Rng + rand_core::CryptoRng + rand::CryptoRng {
    let die = Uniform::from(0..u);
    let val = die.sample(&mut rng);
    Scalar::from(val)
}

pub fn random_scalars_range<R>(mut rng: &mut R, u: u64, n: usize) -> Vec<Scalar> 
    where R: rand_core::RngCore + rand::Rng + rand_core::CryptoRng + rand::CryptoRng {
    
    let mut v = Vec::with_capacity(n);

    for _ in 0..n {
        v.push(random_scalar_range(&mut rng, u));
    }
    v
}


#[derive(Clone)]
pub struct Share {
    pub(crate) share: [Scalar; 2],
}

impl Share {
    pub fn get(&self) -> &[Scalar] {
        self.share.as_slice()
    }
}

/// Checks that the committed degred is low
pub fn low_deg_test(coms: &Vec<G1Projective>, sc: &SharingConfiguration) -> bool {
    // If the degree is n-1, then the check is trivially true
    if sc.t == sc.n {
        return true; 
    }

    let mut rng = thread_rng();
    let batch_dom = BatchEvaluationDomain::new(sc.n);   
    let vf = get_dual_code_word(sc.t - 1, &batch_dom, sc.n, &mut rng);   
    let ip = G1Projective::multi_exp(&coms, vf.as_ref());
    
    ip.eq(&G1Projective::identity())
}

#[allow(unused)]
pub fn get_dual_code_word<R: rand_core::RngCore + rand_core::CryptoRng>(
    deg: usize,
    batch_dom: &BatchEvaluationDomain,
    n: usize,
    mut rng: &mut R,
) -> Vec<Scalar> {
    // The degree-(t-1) polynomial p(X) that shares our secret
    // So, deg = t-1 => t = deg + 1
    // The "dual" polynomial f(X) of degree n - t - 1 = n - (deg + 1) - 1 = n - deg - 2
    let mut f = random_scalars(n - deg - 2, &mut rng);

    // Compute f(\omega^i) for all i's
    let dom = batch_dom.get_subdomain(n);
    fft_assign(&mut f, &dom);
    f.truncate(n);

    // Compute v_i = 1 / \prod_{j \ne i, j \in [0, n-1]} (\omega^i - \omega^j), for all i's
    let v = all_lagrange_denominators(&batch_dom, n);

    // Compute v_i * f(\omega^i), for all i's
    let vf = f
        .iter()
        .zip(v.iter())
        .map(|(v, f)| v.mul(f))
        .collect::<Vec<Scalar>>();

    vf
}