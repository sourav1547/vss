use criterion::{criterion_main, criterion_group};
use criterion::{BenchmarkGroup, BenchmarkId, Criterion, measurement::Measurement, Throughput};
use rand::thread_rng;
use vss::vss::keys::InputSecret;
use vss::pvss::SharingConfiguration;


pub fn vss_ed_group(c: &mut Criterion) {
    let mut group = c.benchmark_group("low-vss-ed");

    let ts = [86, 171, 342];
    let ns= [256, 512, 1024];

    for (&t, &n) in ts.iter().zip(ns.iter()) {
        vss_ed::vss_deal(2*t-1, n, &mut group);
        vss_ed::vss_verify(2*t-1, n, &mut group);
    }

    group.finish();
}


mod vss_ed {
    use super::*;

    use aptos_crypto::Signature;
    use aptos_crypto::ed25519::{Ed25519PrivateKey, Ed25519PublicKey};
    use aptos_crypto::multi_ed25519::MultiEd25519PublicKey;


    use vss::vss::poly_com_ed::{generate_ed_sig_keys, PolyComDealer, PolyComReceiver, Node};
    use vss::vss::public_parameters::PublicParameters;

    #[allow(warnings)]
    pub(crate) fn vss_deal<M: Measurement>(t: usize, n: usize, g: &mut BenchmarkGroup<M>) {
        let mut rng = thread_rng();

        g.throughput(Throughput::Elements(n as u64));

        let pp = PublicParameters::default();
        let sc = SharingConfiguration::new(t, n);
        
        let keys = generate_ed_sig_keys(n);
        let skeys = keys.iter().map(|x| &x.private_key).collect::<Vec<&Ed25519PrivateKey>>();
        let pkeys = keys.iter().map(|x| &x.public_key).collect::<Vec<&Ed25519PublicKey>>();
        

        g.bench_function(BenchmarkId::new(format!("deal-ed-{}", t), n), move |b| {
            b.iter_with_setup(|| {
                let s = InputSecret::new_random(&sc, false, &mut rng);
                let dealer = PolyComDealer::deal(&sc, &pp, &s);
        
                let mut sigs = Vec::with_capacity(n);
                let msg = dealer.msg();

                for i in 0..t {
                    sigs.push(skeys[i].sign_arbitrary_message(msg.as_slice()))
                }
                (s, sigs)
        
            }, |(s, sigs)| {
                let dealer = PolyComDealer::deal(&sc, &pp, &s);
                let msg = dealer.msg();

                for (sig, pk) in sigs.iter().zip(pkeys.iter()) {
                    sig.verify_arbitrary_msg(&msg, pk);
                }
                let mut signers = vec![false; n];
                for i in 0..t {
                    signers[i] = true;
                }
                dealer.get_transcript(&signers, sigs);

            })
        });
    }


    pub(crate) fn vss_verify<M: Measurement>(t: usize, n: usize, g: &mut BenchmarkGroup<M>) {
        let mut rng = thread_rng();

        g.throughput(Throughput::Elements(n as u64));

        let pp = PublicParameters::default();
        let sc = SharingConfiguration::new(t, n);

        let keys = generate_ed_sig_keys(n);
        let skeys = keys.iter().map(|x| &x.private_key).collect::<Vec<&Ed25519PrivateKey>>();
        let pkeys = keys.iter().map(|x| &x.public_key).collect::<Vec<&Ed25519PublicKey>>();

            
        
        g.bench_function(BenchmarkId::new(format!("verify-ed-{}", t), n), move |b| {
            b.iter_with_setup(|| {

                let s = InputSecret::new_random(&sc, true, &mut rng);
                let dealer = PolyComDealer::deal(&sc, &pp, &s);
        
                let mut sigs = Vec::with_capacity(n);
                let msg = dealer.msg();

                for i in 0..t {
                    sigs.push(skeys[i].sign_arbitrary_message(msg.as_slice()))
                }

                let mut signers = vec![false; n];
                for i in 0..t {
                    signers[i] = true;
                }
                let trx = dealer.get_transcript(&signers, sigs);
                let node = Node::new(skeys[0].clone());
                let recv = PolyComReceiver::new(dealer.coms().to_vec(), dealer.shares(0).clone(), node);
                let pkeys2 = pkeys.iter().map(|&x| x.clone()).collect::<Vec<Ed25519PublicKey>>();

                (recv, trx, pkeys2)
            }, |(recv, trx, pkeys2)| {
                recv.sign_verified_deal(&sc, &pp, 0);
             
                let mpk_th = t.try_into().unwrap();
                let mpk = MultiEd25519PublicKey::new(pkeys2, mpk_th).unwrap();
                assert!(recv.verify_transcript(&trx, &sc, &pp, &mpk));
            })
        });
    }
}


criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(10);
    //config = Criterion::default();
    targets = vss_ed_group);
criterion_main!(benches);