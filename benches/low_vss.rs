use criterion::{BenchmarkGroup, BenchmarkId, Criterion, criterion_group, criterion_main, measurement::Measurement, Throughput};
use rand::thread_rng;
use vss::vss::keys::InputSecret;
use vss::pvss::SharingConfiguration;


pub fn vss_group(c: &mut Criterion) {
    let mut group = c.benchmark_group("low-vss");
    
    let ts = [86, 171, 342];
    let ns= [256, 512, 1024];


    for (&t, &n) in ts.iter().zip(ns.iter()) {
        // vss_bls::yurek_deal(t, n, &mut group);
        vss_bls::vss_deal(2*t-1, n, &mut group);
        vss_bls::vss_verify(2*t-1, n, &mut group);
        // vss_bls::yurek_vss_verify(t, n, &mut group);
        // vss_bls::vss_recon(t, n, &mut group);
        // vss_bls::vss_recon(2*t-1, n, &mut group);
    }

    group.finish();
}

mod vss_bls {
    use super::*;

    use aptos_crypto::{SigningKey, Signature};
    use aptos_crypto::bls12381::{PublicKey, PrivateKey};
    use vss::vss::poly_com_bls::{generate_sig_keys, PolyComDealer, PolyComReceiver, Node};
    use vss::vss::public_parameters::PublicParameters;
    use vss::vss::recon::reconstruct;
    use rand::seq::IteratorRandom;


    #[allow(warnings)]
    pub(crate) fn yurek_deal<M: Measurement>(t: usize, n: usize, g: &mut BenchmarkGroup<M>) {
        g.throughput(Throughput::Elements(n as u64));

        let mut rng = thread_rng();
        let pp = PublicParameters::default();
        let sc = SharingConfiguration::new(t, n);
        
        g.bench_function(BenchmarkId::new(format!("yurek-deal-{}", t), n), move |b| {
            b.iter_with_setup(|| {
                let s = InputSecret::new_random(&sc, true, &mut rng);
                s
            }, |s| {
                let dealer = PolyComDealer::deal(&sc, &pp, &s);
            })
        });
    }

    #[allow(warnings)]
    pub(crate) fn vss_recon<M: Measurement>(t: usize, n: usize, g: &mut BenchmarkGroup<M>) {
        g.throughput(Throughput::Elements(n as u64));

        let mut rng = thread_rng();
        let pp = PublicParameters::default();
        let sc = SharingConfiguration::new(t, n);


        g.bench_function(BenchmarkId::new(format!("recon-{}", t), n), move |b| {
            b.iter_with_setup(|| {
                let s = InputSecret::new_random(&sc, true, &mut rng);
                let dealer = PolyComDealer::deal(&sc, &pp, &s);
                
                let mut players : Vec<usize> = (0..n)
                .choose_multiple(&mut rng, t)
                .into_iter().collect::<Vec<usize>>();
                players.sort();

                (s.get_secret_a(), dealer, players)

            }, |(s, dealer, players)| {

                let mut shares = Vec::with_capacity(t);
                for i in 0..t {
                    shares.push(dealer.shares(players[i]));
                }

                let (recon_s, _) = reconstruct(&dealer.coms().to_vec(), &shares, &players, n, &pp);
                assert!(s == recon_s);
            })
        });
    }


    #[allow(warnings)]
    pub(crate) fn yurek_vss_verify<M: Measurement>(t: usize, n: usize, g: &mut BenchmarkGroup<M>) {
        g.throughput(Throughput::Elements(n as u64));

        let mut rng = thread_rng();
        let pp = PublicParameters::default();
        let sc = SharingConfiguration::new(t, n);
        
        g.bench_function(BenchmarkId::new(format!("yurek-verify-{}", t), n), move |b| {
            b.iter_with_setup(|| {
                let s = InputSecret::new_random(&sc, true, &mut rng);
                let dealer = PolyComDealer::deal(&sc, &pp, &s);

                dealer
            }, |dealer| {
                assert!(PolyComReceiver::verify_com(&dealer.coms().to_vec(), &sc)); 
                assert!(PolyComReceiver::verify_eval(&dealer.coms().to_vec(), &pp, 0, &dealer.shares(0)));
            })
        });
    }

    #[allow(warnings)]
    pub(crate) fn vss_deal<M: Measurement>(t: usize, n: usize, g: &mut BenchmarkGroup<M>) {
        let mut rng = thread_rng();

        g.throughput(Throughput::Elements(n as u64));

        let pp = PublicParameters::default();
        let sc = SharingConfiguration::new(t, n);
        
        let keys = generate_sig_keys(n);
        let skeys = keys.iter().map(|x| &x.private_key).collect::<Vec<&PrivateKey>>();
        let pkeys = keys.iter().map(|x| &x.public_key).collect::<Vec<&PublicKey>>();
        

        g.bench_function(BenchmarkId::new(format!("deal-{}", t), n), move |b| {
            b.iter_with_setup(|| {
                let s = InputSecret::new_random(&sc, true, &mut rng);
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

        let keys = generate_sig_keys(n);
        let skeys = keys.iter().map(|x| &x.private_key).collect::<Vec<&PrivateKey>>();
        let pkeys = keys.iter().map(|x| &x.public_key).collect::<Vec<&PublicKey>>();
        
        g.bench_function(BenchmarkId::new(format!("verify-{}", t), n), move |b| {
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
            

                (recv, trx)
            }, |(recv, trx)| {
                recv.sign_verified_deal(&sc, &pp, 0);
                assert!(recv.verify_transcript(&trx, &sc, &pp, &pkeys));
            })
        });
    }


    #[allow(warnings)]
    pub(crate) fn vss_low_deg_test<M: Measurement>(t: usize, n: usize, g: &mut BenchmarkGroup<M>) {
        let mut rng = thread_rng();

        g.throughput(Throughput::Elements(n as u64));

        let pp = PublicParameters::default();
        let sc = SharingConfiguration::new(t, n);

        g.bench_function(BenchmarkId::new(format!("low-deg-test-{}", t), n), move |b| {
            b.iter_with_setup(|| {
                let s = InputSecret::new_random(&sc, true, &mut rng);
                let dealer = PolyComDealer::deal(&sc, &pp, &s);
                let coms = dealer.coms().to_vec();
                coms
            }, |coms| {
                PolyComReceiver::verify_com(&coms, &sc)
            })
        });
    }

}


criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(10);
    //config = Criterion::default();
    targets = vss_group);
criterion_main!(benches);