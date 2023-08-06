use criterion::{BenchmarkGroup, BenchmarkId, Criterion, criterion_group, criterion_main, measurement::Measurement, Throughput};
use rand::thread_rng;
use vss::vss::keys::InputSecret;
use vss::pvss::SharingConfiguration;


pub fn ed_mixed_vss_group(c: &mut Criterion) {
    let mut group = c.benchmark_group("ed-mixed-vss");
    
    let ts = [86, 171, 342];
    let ns= [256, 512, 1024];

    for (&t, &n) in ts.iter().zip(ns.iter()) {
        // vss_mixed_ed::vss_deal(t, n, &mut group);
        vss_mixed_ed::vss_verify(t, n, &mut group);
        // vss_bls::vss_recon(t, n, &mut group);
    }

    group.finish();
}

mod vss_mixed_ed {
    use std::ops::Mul;

    use super::*;

    use aptos_crypto::Signature;
    use aptos_crypto::ed25519::{Ed25519PrivateKey, Ed25519PublicKey};
    use aptos_crypto::multi_ed25519::MultiEd25519PublicKey;


    use vss::vss::poly_com_ed::{generate_ed_sig_keys, Node};
    use vss::vss::poly_com_mixed_ed::{PolyComDealer, PolyComReceiver};
    use vss::vss::public_parameters::PublicParameters;
    use vss::random_scalars;
    use blstrs::G1Projective;
    use group::Group;



    #[allow(warnings)]
    pub(crate) fn vss_deal<M: Measurement>(t: usize, n: usize, g: &mut BenchmarkGroup<M>) {
        let mut rng = thread_rng();

        g.throughput(Throughput::Elements(n as u64));

        let deg = 2*t;
        let pp = PublicParameters::default();
        let sc = SharingConfiguration::new(deg+1, n);
        
        let keys = generate_ed_sig_keys(n);
        let sig_keys = keys.iter().map(|x| &x.private_key).collect::<Vec<&Ed25519PrivateKey>>();
        let ver_keys = keys.iter().map(|x| &x.public_key).collect::<Vec<&Ed25519PublicKey>>();

        let g1 = G1Projective::generator();
        let h = pp.get_bases()[1];
        let dec_keys = random_scalars(n, &mut rng);
        let enc_keys = dec_keys.iter().map(|x| g1.mul(x)).collect::<Vec<_>>();

        g.bench_function(BenchmarkId::new(format!("mixed-deal-ed-{}", t), n), move |b| {
            b.iter_with_setup(|| {
                let s = InputSecret::new_random(&sc, true, &mut rng);
                let dealer = PolyComDealer::deal(&sc, &pp, &s);
        
                let mut sigs = Vec::with_capacity(n);
                let msg = dealer.msg();

                for i in 0..=deg {
                    sigs.push(sig_keys[i].sign_arbitrary_message(msg.as_slice()))
                }
                (s, sigs)
        
            }, |(s, sigs)| {
                let dealer = PolyComDealer::deal(&sc, &pp, &s);
                let msg = dealer.msg();

                for (sig, pk) in sigs.iter().zip(ver_keys.iter()) {
                    sig.verify_arbitrary_msg(&msg, pk);
                }
                let mut signers = vec![false; n];
                for i in 0..=deg {
                    signers[i] = true;
                }

                let tx = dealer.get_transcript(&signers, sigs, &enc_keys, t, deg, &h);
            })
        });
    }


    #[allow(warnings)]
    pub(crate) fn vss_verify<M: Measurement>(t: usize, n: usize, g: &mut BenchmarkGroup<M>) {
        let mut rng = thread_rng();

        g.throughput(Throughput::Elements(n as u64));

        let deg = 2*t;
        let pp = PublicParameters::default();
        let sc = SharingConfiguration::new(deg+1, n);
        
        let keys = generate_ed_sig_keys(n);
        let sig_keys = keys.iter().map(|x| &x.private_key).collect::<Vec<&Ed25519PrivateKey>>();
        let ver_keys = keys.iter().map(|x| &x.public_key).collect::<Vec<&Ed25519PublicKey>>();

        let g1 = G1Projective::generator();
        let h = pp.get_bases()[1];
        let dec_keys = random_scalars(n, &mut rng);
        let enc_keys = dec_keys.iter().map(|x| g1.mul(x)).collect();
        

        g.bench_function(BenchmarkId::new(format!("mixed-ed-verify-{}", t), n), move |b| {
            b.iter_with_setup(|| {
                let s = InputSecret::new_random(&sc, true, &mut rng);
                let dealer = PolyComDealer::deal(&sc, &pp, &s);
        
                let mut sigs = Vec::with_capacity(n);
                let msg = dealer.msg();

                for i in 0..=2*t {
                    sigs.push(sig_keys[i].sign_arbitrary_message(msg.as_slice()))
                }

                let mut signers = vec![false; n];
                for i in 0..=2*t {
                    signers[i] = true;
                }
                let trx = dealer.get_transcript(&signers, sigs, &enc_keys, t, deg, &h);
                let node = Node::new(sig_keys[0].clone());
                let recv = PolyComReceiver::new(dealer.coms().to_vec(), dealer.shares(0).clone(), node);
                let pkeys2 = ver_keys.iter().map(|&x| x.clone()).collect::<Vec<Ed25519PublicKey>>();


                (recv, trx, pkeys2)
        
            }, |(recv, trx, pkeys2)| {
                recv.sign_verified_deal(&sc, &pp, 0);
                
                let mpk_th = t.try_into().unwrap();
                let mpk = MultiEd25519PublicKey::new(pkeys2, mpk_th).unwrap();

                assert!(recv.verify_transcript(&trx, &sc, &pp, &mpk, &enc_keys));
            })
        });
    }


}


criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(10);
    //config = Criterion::default();
    targets = ed_mixed_vss_group);
criterion_main!(benches);