use std::ops::Mul;

use blstrs::{G1Projective, Scalar};
use ff::Field;
use crate::{evaluation_domain::BatchEvaluationDomain, lagrange::lagrange_coefficients_at_zero};
use super::{common::Share, public_parameters::PublicParameters};

pub fn reconstruct(coms: &Vec<G1Projective>, shares: &Vec<Share>, players: &Vec<usize>, n:usize, pp: &PublicParameters) -> (Scalar, Scalar) {
    let batch_dom = BatchEvaluationDomain::new(n);
    let lagr = lagrange_coefficients_at_zero(&batch_dom, players.as_slice());

    let mut s = Scalar::zero();
    let mut r = Scalar::zero();

    let t = shares.len();
    for i in 0..t {
        let com = G1Projective::multi_exp(pp.get_bases(), shares[i].get());
        assert!(coms[players[i]].eq(&com));

        s += lagr[i].mul(shares[i].share[0]);
        r += lagr[i].mul(shares[i].share[1]);
    }

    (s, r)
}