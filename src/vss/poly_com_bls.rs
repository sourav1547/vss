use aptos_bitvec::BitVec;
use aptos_crypto::Uniform;
use aptos_crypto::test_utils::{TEST_SEED, KeyPair};
use blstrs::{G1Projective, Scalar};
use ff::Field;
use rand::rngs::StdRng;
use rand::thread_rng;
use rand_core::SeedableRng;


use crate::vss::common::{random_scalars_range, Share};
use crate::vss::transcript::TranscriptBLS;
use crate::fft::fft;
use crate::vss::public_parameters::PublicParameters;
use crate::vss::keys::InputSecret;
use crate::pvss::SharingConfiguration;

use aptos_crypto::{
    bls12381,
    bls12381::{PrivateKey, PublicKey}, 
    SigningKey, Signature
};

use super::sigs::AggregateSignature;
use super::common::low_deg_test;

pub struct Node {
    pub(crate) sig_key: PrivateKey,
}

impl Node {
    pub fn new(sig_key: PrivateKey) -> Node {
        Self {sig_key}
    }
}

pub struct PolyComDealer {
    /// Total number of nodes
    n : usize,
    /// Pedersen commitment to the polynomial
    pub(crate) coms: Vec<G1Projective>,
    /// Shares of remaining polynomials
    pub(crate) shares: Vec<Share>
}

pub struct PolyComReceiver {
    /// Pedersen commitment to the polynomial
    coms: Vec<G1Projective>,
    /// Shares of remaining polynomials
    share: Share,
    /// PKI information and node index
    node: Node,
}

impl PolyComDealer {

    pub fn deal_coeff(sc: &SharingConfiguration, pp: &PublicParameters, s:&InputSecret) -> Self {
        assert_eq!(s.get_secret_f().len(), sc.t);
        assert_eq!(s.get_secret_r().len(), sc.t);

        let f = s.get_secret_f();
        let r = s.get_secret_r();

        let mut f_evals = fft(f, sc.get_evaluation_domain());
        f_evals.truncate(sc.n);

        let mut r_evals = fft(r, sc.get_evaluation_domain());
        r_evals.truncate(sc.n);

        let mut shares: Vec<Share> = Vec::with_capacity(sc.n);
        for i in 0..sc.n {
            shares.push(Share{share: [f_evals[i], r_evals[i]]});
        }

        let mut coms:Vec<G1Projective> = Vec::with_capacity(sc.t);
        for i in 0..sc.t {
            let scalars = [f[i], r[i]];
            coms.push(G1Projective::multi_exp(pp.get_bases(), scalars.as_slice())); 
        }

        PolyComDealer {
            n: coms.len(),
            coms: coms,
            shares: shares,
        }
    }


    pub fn deal(sc: &SharingConfiguration, pp: &PublicParameters, s:&InputSecret) -> Self {
        assert_eq!(s.get_secret_f().len(), sc.t);
        assert_eq!(s.get_secret_r().len(), sc.t);

        let f = s.get_secret_f();
        let r = s.get_secret_r();

        let mut f_evals = fft(f, sc.get_evaluation_domain());
        f_evals.truncate(sc.n);

        let mut r_evals = fft(r, sc.get_evaluation_domain());
        r_evals.truncate(sc.n);

        let mut shares: Vec<Share> = Vec::with_capacity(sc.n);
        for i in 0..sc.n {
            shares.push(Share{share: [f_evals[i], r_evals[i]]});
        }

        let mut coms:Vec<G1Projective> = Vec::with_capacity(sc.n);
        for i in 0..sc.n {
            let scalars = [f_evals[i], r_evals[i]];
            coms.push(G1Projective::multi_exp(pp.get_bases(), scalars.as_slice())); 
        }

        PolyComDealer {
            n: coms.len(),
            coms: coms,
            shares: shares,
        }
    }

    // Verifies individual signatures
    pub fn verify_sig(&self, pk: &PublicKey, sig: bls12381::Signature) -> bool {
        let msg = bcs::to_bytes(&self.coms).unwrap();
        sig.verify_arbitrary_msg(msg.as_slice(), pk).is_ok()
    }

    // Takes as input a vector of boolean indicating which signers are set
    pub fn aggregate_sig(&self, signers: Vec<bool>, sigs: Vec<bls12381::Signature>) -> AggregateSignature {
        AggregateSignature::new(BitVec::from(signers), Some(bls12381::Signature::aggregate(sigs).unwrap()))
    }

    // This function assumes that all signatures are valid
    pub fn get_transcript(&self, signers: &Vec<bool>, sigs: Vec<bls12381::Signature>) -> TranscriptBLS {
        let agg_sig = self.aggregate_sig(signers.clone(), sigs);
        let missing_count = self.n-agg_sig.get_num_voters();

        let mut shares = Vec::with_capacity(missing_count);
        let mut randomness = Vec::with_capacity(missing_count);

        for (i, &is_set) in signers.iter().enumerate() {
            if !is_set {
                shares.push(self.shares[i].share[0]);
                randomness.push(self.shares[i].share[1]);
            }
        }
        
        TranscriptBLS::new(self.coms.clone(), shares, randomness, agg_sig)
    }

    // TODO: Send hash of the message instead of the entire message.
    pub fn msg(&self) -> Vec<u8> {
        bcs::to_bytes(&self.coms).unwrap()
    }

    pub fn coms(&self) -> &[G1Projective] {
        &self.coms
    }

    pub fn shares(&self, i:usize) -> Share {
        self.shares[i].clone()
    }

}



impl PolyComReceiver {
    pub fn new(coms:Vec<G1Projective>, share:Share, node: Node) -> Self {
        Self {coms, share, node}
    }

    /// Checks whether the commitment is to a low-degree polynomial
    pub fn verify_com(coms:&Vec<G1Projective>, sc: &SharingConfiguration) -> bool {
        low_deg_test(coms, sc)
    }

    /// Checks correct evaluation at a point among `[1,...n]`
    pub fn verify_eval(coms:&Vec<G1Projective>, pp: &PublicParameters, i:usize, share: &Share) -> bool {
        let com = G1Projective::multi_exp(pp.get_bases(), share.get());
        coms[i].eq(&com)
    }

    pub fn sign_verified_deal(&self, sc:&SharingConfiguration, pp: &PublicParameters, i:usize) -> Option<bls12381::Signature> {
        // Return signature the dealing is valid
        if PolyComReceiver::verify_com(&self.coms, sc) & PolyComReceiver::verify_eval(&self.coms, pp, i, &self.share) {
            let msg =  bcs::to_bytes(&self.coms).unwrap();
            return Some(self.node.sig_key.sign_arbitrary_message(msg.as_slice()));
        }
        None
    }
    
    pub fn verify_transcript(&self, t: &TranscriptBLS, sc: &SharingConfiguration, pp: &PublicParameters, pks: &Vec<&PublicKey>) -> bool {
        let num_signed = t.agg_sig().get_num_voters();
        let n = t.coms().len();
        let missing_ct = n-num_signed;
        
        // Checking low-degree of the committed polynomial
        assert!(PolyComReceiver::verify_com(t.coms(), sc)); 
        assert!(t.shares().len() == t.randomness().len());
        assert!(t.shares().len() == missing_ct);

        // TODO: To interpolate the commitments at 0 in the exponent and check that it matches with the t.g_at_0

        // Checking correctness of aggregate signature
        let msg = bcs::to_bytes(&self.coms).unwrap();
        assert!(t.agg_sig().verify(msg.as_slice(), pks));

        // Encryption keys of nodes whose shares are not opened
        // FIXME: Currently this is just a placeholder
        let mut missing_coms = Vec::with_capacity(missing_ct);
        let mut rng = thread_rng();
        let lambdas = random_scalars_range(&mut rng, u64::MAX, missing_ct);

        // Checking the correctness of the revealed shares and randomness 
        let mut idx = 0;
        let mut s = Scalar::zero();
        let mut r = Scalar::zero();
        for pos in 0..n {
            if !t.agg_sig().get_signers_bitvec().is_set(pos as u16) {
                s += lambdas[idx]*t.shares()[idx];
                r += lambdas[idx]*t.randomness()[idx];
                
                idx +=1;
                missing_coms.push(self.coms[pos]);
            }
        }

        let com_pos = G1Projective::multi_exp(pp.get_bases(), [s, r].as_slice());
        let com = G1Projective::multi_exp(&missing_coms, &lambdas);
        
        com_pos == com
    }

    pub fn get_share(&self) -> Share {
        self.share.clone()
    }
}


// Helper function to generate N bls12381 private keys.
pub fn generate_sig_keys(n: usize) -> Vec<KeyPair<PrivateKey, PublicKey>> {
    let mut rng = StdRng::from_seed(TEST_SEED);
    (0..n)
        .map(|_| KeyPair::<PrivateKey, PublicKey>::generate(&mut rng))
        .collect()
}


#[allow(unused)]
mod test {
    use super::*;
    use group::Group;
    use std::ops::Mul;
    use rand::thread_rng;
    use rand::seq::IteratorRandom;
    use crate::vss::recon::reconstruct;
    use crate::random_scalars;

    #[test]
    pub fn test_vss(){
        let mut rng = thread_rng();

        // for t in [1, 2, 3, 4, 5, 6, 7, 8] {
        for th in [2, 4, 6, 8] {
            let n = 3*(th-1) + 1;
            // let n = th+1;
            // let th = 2*th -1;
            let pp = PublicParameters::default();
            let sc = SharingConfiguration::new(th, n);
            let keys = generate_sig_keys(n);
            let ver_keys = keys.iter().map(|x| &x.public_key).collect::<Vec<&PublicKey>>();

            let s = InputSecret::new_random(&sc, true, &mut rng);
            let dealer = PolyComDealer::deal(&sc, &pp, &s);

            let mut sigs = Vec::with_capacity(n);
            let mut recvs = vec![];

            for i in 0..n {
                let node = Node::new(keys[i].private_key.clone());

                let recv = PolyComReceiver::new(dealer.coms.clone(), dealer.shares[i].clone(), node);
                let sig = recv.sign_verified_deal(&sc, &pp, i);
                assert!(sig.is_some());
                
                sigs.push(sig.unwrap());
                recvs.push(recv);
            }

            let mut signers = vec![false; n];
            let mut valid_sigs: Vec<bls12381::Signature> = Vec::with_capacity(th);

            let msg = dealer.msg();
            let mut num_valid =0;
            for i in 0..n {
                if sigs[i].verify_arbitrary_msg(msg.as_slice(), ver_keys[i]).is_ok() {
                    signers[i] = true;
                    valid_sigs.push(sigs[i].clone());
                    num_valid = num_valid+1;
                    if num_valid == th {
                        break
                    }
                }
            }
            let t = dealer.get_transcript(&signers, valid_sigs);
            
            assert!(t.agg_sig().verify(&msg, &ver_keys));

            // TODO: I probably have to check marshalling as well.
            for recv in recvs.iter() {
                assert!(recv.verify_transcript(&t, &sc, &pp, &ver_keys));
            }

            let mut shares = Vec::with_capacity(th);
            
            let mut players : Vec<usize> = (0..n)
            .choose_multiple(&mut rng, th)
            .into_iter().collect::<Vec<usize>>();
            players.sort();

            for i in 0..th {
                shares.push(dealer.shares(players[i]));
                // shares.push(recvs[players[i]].get_share());
            }

            let (recon_s, _) = reconstruct(&dealer.coms, &shares, &players, n, &pp);
            assert!(s.get_secret_a() == recon_s);
        }
    }

}