use aptos_bitvec::BitVec;
use aptos_crypto::{Uniform, Signature};
use aptos_crypto::ed25519::{Ed25519PrivateKey, Ed25519PublicKey, Ed25519Signature};
use aptos_crypto::multi_ed25519::{MultiEd25519PublicKey, MultiEd25519Signature};
use aptos_crypto::test_utils::{TEST_SEED, KeyPair};
use blstrs::{G1Projective, Scalar};
use ff::Field;
use rand::rngs::StdRng;
use rand::thread_rng;
use rand_core::SeedableRng;

use crate::vss::common::random_scalars_range;
use crate::vss::transcript::TranscriptEd;
use crate::fft::fft;
use crate::vss::public_parameters::PublicParameters;
use crate::vss::keys::InputSecret;
use crate::pvss::SharingConfiguration;

use super::sigs::EdSignature;
use super::common::{Share, low_deg_test};

pub struct Node {
    pub(crate) sig_key: Ed25519PrivateKey,
}

impl Node {
    pub fn new(sig_key:Ed25519PrivateKey) -> Self {
        Self {sig_key}
    }
}

pub struct PolyComDealer {
    /// Total number of nodes
    n : usize,
    /// Pedersen commitment to the polynomial
    coms: Vec<G1Projective>,
    /// Shares of remaining polynomials
    shares: Vec<Share>
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
    pub fn deal(sc: &SharingConfiguration, pp: &PublicParameters, s: &InputSecret) -> Self {
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
    pub fn verify_sig(&self, pk: &Ed25519PublicKey, sig: Ed25519Signature) -> bool {
        let msg = bcs::to_bytes(&self.coms).unwrap();
        sig.verify_arbitrary_msg(msg.as_slice(), pk).is_ok()
    }

    // Takes as input a vector of boolean indicating which signers are set
    pub fn aggregate_sig(&self, signers: Vec<bool>, sigs: Vec<Ed25519Signature>) -> EdSignature {
        // AggregateSignature::new(BitVec::from(signers), Some(bls12381::Signature::aggregate(sigs).unwrap()))
        let mut indices: Vec<usize> = Vec::with_capacity(sigs.len());
        for i in 0..signers.len() {
            if signers[i] {
                indices.push(i);
            }
        }

        let new_sigs = sigs.iter().zip(indices.iter()).map(|(s, &i)| (s.clone(),i)).collect::<Vec<(Ed25519Signature,usize)>>();
        let mt_sig = MultiEd25519Signature::new(new_sigs);
        EdSignature::new(BitVec::from(signers), Some(mt_sig.unwrap()))
    }

    // This function outputs the Mixed-VSS transcript. 
    // This function assumes that all signatures are valid
    pub fn get_transcript(&self, signers: &Vec<bool>, sigs: Vec<Ed25519Signature>) -> TranscriptEd {
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

        TranscriptEd::new(self.coms.clone(), shares, randomness, agg_sig)

    }

    // TODO: Send hash of the message instead of the entire message.
    pub fn msg(&self) -> Vec<u8> {
        bcs::to_bytes(&self.coms).unwrap()
    }

    pub fn coms(&self) -> &[G1Projective] {
        &self.coms
    }

    pub fn shares(&self, i:usize) -> &Share {
        &self.shares[i]
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

    pub fn verify_eval(coms:&Vec<G1Projective>, pp: &PublicParameters, i:usize, share: &Share) -> bool {
        let com = G1Projective::multi_exp(pp.get_bases(), share.get());
        coms[i].eq(&com)
    }

    pub fn sign_verified_deal(&self, sc:&SharingConfiguration, pp: &PublicParameters, i:usize) -> Option<Ed25519Signature> {
        // Return signature the dealing is valid
        if PolyComReceiver::verify_com(&self.coms, sc) & PolyComReceiver::verify_eval(&self.coms, pp, i, &self.share) {
            let msg =  bcs::to_bytes(&self.coms).unwrap();
            return Some(self.node.sig_key.sign_arbitrary_message(msg.as_slice()));
        }
        None
    }
    
    pub fn verify_transcript(&self, t: &TranscriptEd, sc: &SharingConfiguration, pp: &PublicParameters, pk: &MultiEd25519PublicKey) -> bool {
        let num_signed = t.agg_sig().get_num_voters();
        let n = t.coms().len();
        let missing_ct = n-num_signed;
        
        // Checking low-degree of the committed polynomial
        assert!(PolyComReceiver::verify_com(t.coms(), sc)); 
        assert!(t.shares().len() == t.randomness().len());
        assert!(t.shares().len() == missing_ct);

        // Checking correctness of aggregate signature
        let msg = bcs::to_bytes(&self.coms).unwrap();
        assert!(t.agg_sig().verify(msg.as_slice(), &pk));

        let mut missing_coms = Vec::with_capacity(t.shares().len());

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
}


pub fn generate_ed_sig_keys(n: usize) -> Vec<KeyPair<Ed25519PrivateKey, Ed25519PublicKey>> {
    let mut rng = StdRng::from_seed(TEST_SEED);
    (0..n)
        .map(|_| KeyPair::<Ed25519PrivateKey, Ed25519PublicKey>::generate(&mut rng))
        .collect()
}

#[allow(unused)]
mod test {
    use super::*;
    use group::Group;
    use std::ops::Mul;
    use rand::thread_rng;
    use rand::seq::IteratorRandom;
    use crate::random_scalars;

    #[test]
    pub fn test_vss(){
        let mut rng = thread_rng();

        for th in [2, 4, 8] {
        // for th in [86] {
            let n = 3*th + 1;
            // let n = th+1;

            let pp = PublicParameters::default();
            let sc = SharingConfiguration::new(th, n);
            let keys = generate_ed_sig_keys(n);
            let ver_keys = keys.iter().map(|x| &x.public_key).collect::<Vec<&Ed25519PublicKey>>();
            
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
            let mut valid_sigs = Vec::with_capacity(th);

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
            
            let pub_keys_2 = ver_keys.iter().map(|&x| x.clone()).collect::<Vec<Ed25519PublicKey>>();
            let mpk_th = th.try_into().unwrap();

            let mpk = MultiEd25519PublicKey::new(pub_keys_2, mpk_th).unwrap();
            assert!(t.agg_sig().verify(&msg, &mpk));

            // TODO: I probably have to check marshalling as well.
            for recv in recvs.iter() {
                assert!(recv.verify_transcript(&t, &sc, &pp, &mpk));
            }

        }
    }
}