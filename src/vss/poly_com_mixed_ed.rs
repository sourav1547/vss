use std::cmp;

use aptos_bitvec::BitVec;
use aptos_crypto::Signature;
use aptos_crypto::ed25519::{Ed25519Signature, Ed25519PublicKey};
use aptos_crypto::multi_ed25519::{MultiEd25519Signature, MultiEd25519PublicKey};
use blstrs::{G1Projective, Scalar};


use crate::vss::ni_vss::dealing::verify_dealing;
use crate::vss::ni_vss::encryption::dec_chunks;
use crate::vss::transcript::TranscriptMixedEd;
use crate::fft::fft;
use crate::vss::public_parameters::PublicParameters;
use crate::vss::keys::InputSecret;
use crate::pvss::SharingConfiguration;


use super::poly_com_ed::Node;
use super::sigs::EdSignature;
use super::common::{Share, low_deg_test};
use super::ni_vss::dealing::create_dealing;
use super::ni_vss::encryption::CiphertextChunks;
use super::ni_vss::nizk_chunking::ProofChunking;
use super::ni_vss::nizk_sharing::ProofSharing;


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
    pub fn groth_deal(sc: &SharingConfiguration, pp: &PublicParameters, enc_keys: &[G1Projective], s: &InputSecret)
    -> (Vec<G1Projective>, CiphertextChunks, G1Projective, Vec<G1Projective>, ProofChunking, ProofSharing)
    {
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

        let h = pp.get_bases()[1];
        let (ctxt, r_bb, enc_rr, enc_pf, sh_pf) = create_dealing(&h, &coms, enc_keys, &f_evals, &r_evals);
        
        (coms, ctxt, r_bb, enc_rr, enc_pf, sh_pf)
        
    }


    pub fn groth_verify(
        h: &G1Projective, 
        enc_coms: &[G1Projective], 
        enc_keys: &[G1Projective],
        ciphertext: &CiphertextChunks,
        chunk_pf: &ProofChunking,
        r_bb: &G1Projective, 
        enc_rr: &Vec<G1Projective>,
        share_pf: &ProofSharing,
        dec_key: Scalar,
        idx: usize,
    ) {
        assert!(verify_dealing(h, enc_coms, enc_keys, ciphertext, chunk_pf, r_bb, enc_rr, share_pf));
        dec_chunks(ciphertext, dec_key, idx);
    }

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
        
        let mut indices: Vec<usize> = Vec::with_capacity(sigs.len());
        for i in 0..signers.len() {
            if signers[i] {
                indices.push(i);
            }
        }

        let new_sigs = sigs.iter().zip(indices.iter()).map(|(s, &i)| (s.clone(),i)).collect::<Vec<(Ed25519Signature, usize)>>();
        let mt_sig = MultiEd25519Signature::new(new_sigs);
        EdSignature::new(BitVec::from(signers), Some(mt_sig.unwrap()))
    }

    // This function outputs the Mixed-VSS transcript. 
    // This function assumes that all signatures are valid
    pub fn get_transcript(&self, signers: &Vec<bool>, sigs: Vec<Ed25519Signature>, public_keys: &Vec<G1Projective>, th: usize, deg: usize, h: &G1Projective) -> TranscriptMixedEd {
        assert!(sigs.len() >= 2*th+1);
        let max_reveal = 2*th - deg;
        
        let missing_count = self.n - sigs.len();
        let reveal_count = cmp::min(missing_count, max_reveal);
        let enc_count = cmp::max(0, missing_count - reveal_count);

        let mut reveal_shares = Vec::with_capacity(reveal_count);
        let mut reveal_randomness = Vec::with_capacity(reveal_count);
        
        let mut enc_shares = Vec::with_capacity(enc_count);
        let mut enc_randomness = Vec::with_capacity(enc_count);
        let mut enc_commits = Vec::with_capacity(enc_count);
        let mut enc_pks = Vec::with_capacity(enc_count);


        let mut count = 0;
        for (i, &is_set) in signers.iter().enumerate() {
            if !is_set {
                if count < reveal_count {
                    reveal_shares.push(self.shares[i].share[0]);
                    reveal_randomness.push(self.shares[i].share[1]);
                } else {
                    enc_shares.push(self.shares[i].share[0]);
                    enc_randomness.push(self.shares[i].share[1]);
                    enc_commits.push(self.coms[i]);
                    enc_pks.push(public_keys[i]);
                }
                count +=1;
            }
        }
        
        let agg_sig = self.aggregate_sig(signers.clone(), sigs);
        // let (ciphertext, proof) = create_dealing(&enc_shares, &enc_randomness, &enc_pks);

        let (ciphertext, r_bb, enc_rr, chunk_pf, share_pf) = create_dealing(h, &enc_commits, &enc_pks, &enc_shares, &enc_randomness);
        
        TranscriptMixedEd::new(self.coms.clone(), reveal_shares, reveal_randomness, agg_sig, ciphertext, chunk_pf, r_bb, enc_rr, share_pf)
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
    pub fn verify_eval(&self, pp: &PublicParameters, i:usize, share: &Share) -> bool {
        // let gs = [self.g, self.h].as_slice();
        let com = G1Projective::multi_exp(pp.get_bases(), share.get());
        self.coms[i].eq(&com)
    }

    pub fn sign_verified_deal(&self, sc:&SharingConfiguration, pp: &PublicParameters, i:usize) -> Option<Ed25519Signature> {
        // Return signature the dealing is valid
        if PolyComReceiver::verify_com(&self.coms, sc) & self.verify_eval(pp, i, &self.share) {
            let msg =  bcs::to_bytes(&self.coms).unwrap();
            return Some(self.node.sig_key.sign_arbitrary_message(msg.as_slice()));
        }
        None
    }
    
    pub fn verify_transcript(&self, t: &TranscriptMixedEd, sc: &SharingConfiguration, pp: &PublicParameters, pk: &MultiEd25519PublicKey, pub_keys: &Vec<G1Projective>) -> bool {
        let n = t.coms().len();
        let num_signed = t.agg_sig().get_num_voters();
        let missing_count = n-num_signed;
        let reveal_count = t.reveal_count();
        let enc_count = missing_count-reveal_count;
        
        // Checking low-degree of the committed polynomial
        assert!(PolyComReceiver::verify_com(t.coms(), sc)); 

        // Checking lengths of the openning vectors
        assert!(t.randomness().len() == reveal_count);
        assert!(t.enc_rr().len() == enc_count);

        // Checking correctness of aggregate signature
        let msg = bcs::to_bytes(&self.coms).unwrap();
        assert!(t.agg_sig().verify(msg.as_slice(), &pk));

        // Encryption keys of nodes whose shares are not opened
        // FIXME: Currently this is just a placeholder
        let mut enc_coms = Vec::with_capacity(enc_count);
        let mut enc_keys = Vec::with_capacity(enc_count);
        

        // Checking the correctness of the revealed shares and randomness 
        let mut idx = 0;
        for pos in 0..n {
            if !t.agg_sig().get_signers_bitvec().is_set(pos as u16) {
                
                if idx < reveal_count {
                    let s = t.shares()[idx];
                    let r = t.randomness()[idx];

                    let com_pos = G1Projective::multi_exp(pp.get_bases(), [s, r].as_slice());
                    assert!(com_pos == self.coms[pos]);
                } else {
                    enc_coms.push(self.coms[pos].clone());
                    enc_keys.push(pub_keys[pos].clone());
                }
                idx +=1;
            }
        }

        // FIXME: As of now this assert will fail if the randomness vector is non-zero
        let h = pp.get_bases()[1];
        assert!(verify_dealing(&h, &enc_coms, &enc_keys, &t.ciphertext, &t.chunk_pf, &t.r_bb, &t.enc_rr, &t.share_pf));

        true
    }

    pub fn get_share(&self) -> Share {
        self.share.clone()
    }
}

#[allow(unused)]
mod test {

    use super::*;
    use group::Group;
    use std::ops::Mul;
    use rand::thread_rng;
    use rand::seq::IteratorRandom;
    use crate::vss::poly_com_ed::generate_ed_sig_keys;
    use crate::vss::recon::reconstruct;
    use crate::random_scalars;


    #[test]
    pub fn test_vss(){
        
        let mut rng = thread_rng();

        for th in [4, 8] {
            let n = 3*th + 1; 

            let deg = 2*th;
            
            let pp = PublicParameters::default();
            let sc = SharingConfiguration::new(deg+1, n);
            
            let keys = generate_ed_sig_keys(n);
            let ver_keys = keys.iter().map(|x| &x.public_key).collect::<Vec<&Ed25519PublicKey>>();
            

            let g = G1Projective::generator();
            let dec_keys = random_scalars(n, &mut rng);
            let enc_keys = dec_keys.iter().map(|x| g.mul(x)).collect::<Vec<_>>();

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
                    if num_valid == 2*th+1 {
                        break
                    }
                }
            }
            let h = pp.get_bases()[1];
            let t = dealer.get_transcript(&signers, valid_sigs, &enc_keys, th, deg, &h);

            let pub_keys_2 = ver_keys.iter().map(|&x| x.clone()).collect::<Vec<Ed25519PublicKey>>();
            let mpk_th = th.try_into().unwrap();

            let mpk = MultiEd25519PublicKey::new(pub_keys_2, mpk_th).unwrap();
            assert!(t.agg_sig().verify(&msg, &mpk));


            // TODO: I probably have to check marshalling as well.
            for recv in recvs.iter() {
                assert!(recv.verify_transcript(&t, &sc, &pp, &mpk, &enc_keys));
            }

            let mut shares = Vec::with_capacity(th);
            let mut players : Vec<usize> = (0..n)
            .choose_multiple(&mut rng, deg+1)
            .into_iter().collect::<Vec<usize>>();
            players.sort();

            for i in 0..=deg {
                shares.push(recvs[players[i]].get_share());
            }

            let (recon_s, _) = reconstruct(&dealer.coms, &shares, &players, n, &pp);
            assert!(s.get_secret_a() == recon_s);
        }
    }

}