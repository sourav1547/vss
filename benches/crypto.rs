use std::ops::Mul;
use blstrs::{G1Projective, G2Projective, Gt};
use criterion::{BenchmarkGroup, BenchmarkId, Criterion, criterion_group, criterion_main, measurement::Measurement, Throughput};
use rand::{Rng, thread_rng};
use rand_core::RngCore;

use vss::{hash_to_scalar, polynomials, random_g1_point, random_g1_points, random_g2_point, random_gt_point, random_gt_points, random_scalar, random_scalars};
use vss::evaluation_domain::{BatchEvaluationDomain, EvaluationDomain};
use vss::fft::fft_assign;
use vss::tests::{FFT_SIZES, LARGE_SIZES, OUR_THRESHOLD, SMALL_SIZES};

pub fn crypto_group(c: &mut Criterion) {
    let mut group = c.benchmark_group("crypto");

    batch_evaluation_domain_new(OUR_THRESHOLD, &mut group);
    fft_assign_bench(OUR_THRESHOLD, &mut group);

    gt_multiexp_naive(OUR_THRESHOLD, &mut group);
    g1_multiexp(OUR_THRESHOLD, &mut group);

    accumulator_poly(OUR_THRESHOLD, &mut group);
    accumulator_poly_slow(OUR_THRESHOLD, &mut group);
    random_scalars_and_points_benches(&mut group);

    // Tried to find some optimal parameters but no luck
    for naive in [16, 32, 64, 128, 256, 512] {
        for fft in [16, 32, 64, 128, 256, 512] {
            if fft / naive > 2 {
                continue
            }
            accumulator_poly_scheduled(OUR_THRESHOLD, naive, fft, &mut group);
        }
    }

    for n in LARGE_SIZES {
        g1_multiexp(n, &mut group);
        fft_assign_bench(n, &mut group);
    }

    for n in SMALL_SIZES {
        hash_to_scalar_bench(n, &mut group);
        hash_to_g1_bench(n, &mut group);
        hash_to_g2_bench(n, &mut group);
        evaluation_domain_new(n, &mut group);
        batch_evaluation_domain_get_subdomain(n, &mut group);
        poly_mul_slow(n, &mut group);
        poly_mul_less_slow(n, &mut group);
    }

    for n in FFT_SIZES {
        poly_mul_fft(n, &mut group);
        poly_mul_fft_with_dom(n, &mut group);
    }

    group.finish();
}


fn random_scalars_and_points_benches<M: Measurement>(g: &mut BenchmarkGroup<M>) {
    let mut rng = thread_rng();

    g.bench_function("random_scalar", move |b| {
        b.iter(|| {
            random_scalar(&mut rng)
        })
    });

    g.bench_function("random_g1", move |b| {
        b.iter(|| {
            random_g1_point(&mut rng)
        })
    });

    g.bench_function("random_g2", move |b| {
        b.iter(|| {
            random_g2_point(&mut rng)
        })
    });

    g.bench_function("random_gt", move |b| {
        b.iter(|| {
            random_gt_point(&mut rng)
        })
    });
}

fn hash_to_scalar_bench<M: Measurement>(n: usize, g: &mut BenchmarkGroup<M>) {
    let mut rng = thread_rng();

    g.bench_function(BenchmarkId::new("hash_to_scalar", n), move |b| {
        b.iter_with_setup(|| {
            let mut bytes = Vec::with_capacity(n);
            bytes.resize(n, 0u8);

            rng.fill(bytes.as_mut_slice());

            bytes
        }, |bytes| {
            hash_to_scalar(bytes.as_slice(), b"criterion benchmark domain separation tag")
        })
    });
}

fn hash_to_g1_bench<M: Measurement>(n: usize, g: &mut BenchmarkGroup<M>) {
    let mut rng = thread_rng();

    g.bench_function(BenchmarkId::new("hash_to_g1", n), move |b| {
        b.iter_with_setup(|| {
            let mut bytes = Vec::with_capacity(n);
            bytes.resize(n, 0u8);

            rng.fill(bytes.as_mut_slice());

            bytes
        }, |bytes| {
            G1Projective::hash_to_curve(bytes.as_slice(),b"criterion benchmark domain separation tag", b"criterion benchmark augmented data");
        })
    });
}

fn hash_to_g2_bench<M: Measurement>(n: usize, g: &mut BenchmarkGroup<M>) {
    let mut rng = thread_rng();

    g.bench_function(BenchmarkId::new("hash_to_g2", n), move |b| {
        b.iter_with_setup(|| {
            let mut bytes = Vec::with_capacity(n);
            bytes.resize(n, 0u8);

            rng.fill(bytes.as_mut_slice());

            bytes
        }, |bytes| {
            G2Projective::hash_to_curve(bytes.as_slice(),b"criterion benchmark domain separation tag", b"criterion benchmark augmented data");
        })
    });
}

fn batch_evaluation_domain_new<M: Measurement>(n: usize, g: &mut BenchmarkGroup<M>) {
    g.bench_function( BenchmarkId::new("batch_evaluation_domain::new", n), move |b| {
        b.iter(|| {
            BatchEvaluationDomain::new(n)
        })
    });
}

#[allow(non_snake_case)]
fn batch_evaluation_domain_get_subdomain<M: Measurement>(N: usize, g: &mut BenchmarkGroup<M>) {
    let mut rng = thread_rng();

    g.bench_function(BenchmarkId::new("batch_evaluation_domain::get_subdomain", N), move |b| {
        b.iter_with_setup(|| {
            let batch_dom = BatchEvaluationDomain::new(N);
            let k = rng.next_u64() as usize % batch_dom.N() + 1;
            (batch_dom, k)
        }, |(batch_dom, k)| {
            batch_dom.get_subdomain(k);
        })
    });
}

fn evaluation_domain_new<M: Measurement>(n: usize, g: &mut BenchmarkGroup<M>) {
    g.bench_function( BenchmarkId::new("evaluation_domain::new", n), move |b| {
        b.iter(|| {
            EvaluationDomain::new(n)
        })
    });
}

fn fft_assign_bench<M: Measurement>(n: usize, g: &mut BenchmarkGroup<M>) {
    let mut rng = thread_rng();

    g.throughput(Throughput::Elements(n as u64));

    g.bench_function(BenchmarkId::new("fft_assign", n), move |b| {
        b.iter_with_setup(|| {
            let poly = random_scalars(n, &mut rng);
            let dom = EvaluationDomain::new(n).unwrap();
            (poly, dom)
        }, |(mut poly, dom)| {
            fft_assign(&mut poly, &dom);
        })
    });
}

fn poly_mul_fft<M: Measurement>(n: usize, g: &mut BenchmarkGroup<M>) {
    let mut rng = thread_rng();

    g.throughput(Throughput::Elements(n as u64));

    g.bench_function(BenchmarkId::new("poly_mul_assign_fft", n), move |b| {
        b.iter_with_setup(|| {
            let f = random_scalars(n, &mut rng);
            let g = random_scalars(n, &mut rng);

            (f,g)
        }, |(mut f, mut g)| {
            polynomials::poly_mul_assign_fft(&mut f, &mut g);
        })
    });
}

fn poly_mul_fft_with_dom<M: Measurement>(n: usize, g: &mut BenchmarkGroup<M>) {
    let mut rng = thread_rng();

    g.throughput(Throughput::Elements(n as u64));

    g.bench_function(BenchmarkId::new("poly_mul_assign_fft_with_dom", n), move |b| {
        b.iter_with_setup(|| {
            let f = random_scalars(n, &mut rng);
            let g = random_scalars(n, &mut rng);

            (f,g, EvaluationDomain::new(2*n-1).unwrap())
        }, |(mut f, mut g, dom)| {
            polynomials::poly_mul_assign_fft_with_dom(&mut f, &mut g, &dom);
        })
    });
}

fn poly_mul_slow<M: Measurement>(n: usize, g: &mut BenchmarkGroup<M>) {
    let mut rng = thread_rng();

    g.throughput(Throughput::Elements(n as u64));

    g.bench_function(BenchmarkId::new("poly_mul_slow", n), move |b| {
        b.iter_with_setup(|| {
            let f = random_scalars(n, &mut rng);
            let g = random_scalars(n, &mut rng);

            (f,g)
        }, |(f, g)| {
            polynomials::poly_mul_slow(&f, &g);
        })
    });
}

fn poly_mul_less_slow<M: Measurement>(n: usize, g: &mut BenchmarkGroup<M>) {
    assert_eq!((n & (n - 1)), 0);   // should be a power of two

    let mut rng = thread_rng();

    g.throughput(Throughput::Elements(n as u64));

    g.bench_function(BenchmarkId::new("poly_mul_less_slow", n), move |b| {
        b.iter_with_setup(|| {
            let f = random_scalars(n, &mut rng);
            let g = random_scalars(n, &mut rng);

            (f,g)
        }, |(f, g)| {
            polynomials::poly_mul_less_slow(&f, &g);
        })
    });
}

#[allow(non_snake_case)]
fn accumulator_poly<M: Measurement>(n: usize, g: &mut BenchmarkGroup<M>) {
    let mut rng = thread_rng();
    let FFT_THRESH = 128; // 256 seems to produce the same result

    g.throughput(Throughput::Elements(n as u64));

    g.bench_function(BenchmarkId::new("accumulator_poly", n), move |b| {
        b.iter_with_setup(|| {
            (random_scalars(n, &mut rng), BatchEvaluationDomain::new(n))
        }, |(set, batch_dom)| {
            polynomials::accumulator_poly(set.as_slice(), &batch_dom, FFT_THRESH);
        })
    });
}

fn accumulator_poly_slow<M: Measurement>(n: usize, g: &mut BenchmarkGroup<M>) {
    let mut rng = thread_rng();

    g.throughput(Throughput::Elements(n as u64));

    g.bench_function(BenchmarkId::new("accumulator_poly_slow", n), move |b| {
        b.iter_with_setup(|| {
            random_scalars(n, &mut rng)
        }, |set | {
            polynomials::accumulator_poly_slow(set.as_slice());
        })
    });
}

#[allow(non_snake_case)]
fn accumulator_poly_scheduled<M: Measurement>(n: usize, naive_thresh: usize, fft_thresh: usize, g: &mut BenchmarkGroup<M>) {
    let mut rng = thread_rng();

    // 16 FFT, 32 naive -> 24.3ms
    // 64 FFT, 128 naive -> 17ms
    // 128 FFT, 256 naive -> 15.4 ms
    // 256 FFT, 512 naive -> 14.8 ms
    // 256 FFT, 128 naive -> 14.1 ms

    g.throughput(Throughput::Elements(n as u64));

    let name = format!("accumulator_poly_scheduled/{}/naive={}/fft={}", n, naive_thresh, fft_thresh);
    g.bench_function(name, move |b| {
        b.iter_with_setup(|| {
            (random_scalars(n, &mut rng), BatchEvaluationDomain::new(n))
        }, |(set, batch_dom)| {
            polynomials::accumulator_poly_scheduled(set.as_slice(), &batch_dom, naive_thresh, fft_thresh);
        })
    });
}

fn g1_multiexp<M: Measurement>(n: usize, g: &mut BenchmarkGroup<M>) {
    let mut rng = thread_rng();

    g.throughput(Throughput::Elements(n as u64));

    g.bench_function(BenchmarkId::new("g1_multiexp", n), move |b| {
        b.iter_with_setup(|| {
            let points = random_g1_points(n, &mut rng);

            let scalars = random_scalars(n, &mut rng);

            (points, scalars)
        }, |(points, scalars)| {
            G1Projective::multi_exp(points.as_slice(), scalars.as_ref());
        })
    });
}


fn gt_multiexp_naive<M: Measurement>(n: usize, g: &mut BenchmarkGroup<M>) {
    let mut rng = thread_rng();

    g.throughput(Throughput::Elements(n as u64));

    g.bench_function(BenchmarkId::new("gt_multiexp_naive", n), move |b| {
        b.iter_with_setup(|| {
            let points = random_gt_points(n, &mut rng);

            let scalars = random_scalars(n, &mut rng);

            (points, scalars)
        }, |(points, scalars)| {
            points.into_iter().zip(scalars.into_iter()).map(|(p, s)| p.mul(s)).sum::<Gt>()
        })
    });
}

criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(10);
    //config = Criterion::default();
    targets = crypto_group);
criterion_main!(benches);
