use criterion::{BenchmarkGroup, BenchmarkId, Criterion, criterion_group, criterion_main, measurement::Measurement, Throughput};
use more_asserts::{assert_ge, assert_le};
use rand::seq::IteratorRandom;
use rand::thread_rng;
use vss::evaluation_domain::BatchEvaluationDomain;
use vss::lagrange::lagrange_coefficients_at_zero;
use vss::tests::{OUR_N, OUR_THRESHOLD};

pub fn lagrange_group(c: &mut Criterion) {
    let mut group = c.benchmark_group("lagrange");

    lagrange_tcz20(OUR_THRESHOLD, OUR_N, &mut group);

    group.finish();
}

#[allow(non_snake_case)]
pub fn lagrange_tcz20<M: Measurement>(thresh: usize, n: usize, g: &mut BenchmarkGroup<M>) {
    assert_ge!(thresh, 1);
    assert_le!(thresh, n);
    let mut rng = thread_rng();

    g.throughput(Throughput::Elements(n as u64));

    g.bench_function(BenchmarkId::new(format!("tcz20-thresh={thresh}"), n), move |b| {
        b.iter_with_setup(|| {
            let players : Vec<usize> = (0..n)
                .choose_multiple(&mut rng, thresh)
                .into_iter().collect::<Vec<usize>>();

            let batch_dom = BatchEvaluationDomain::new(n);

            (players, batch_dom)
        }, |(players, batch_dom)| {
            lagrange_coefficients_at_zero(&batch_dom, players.as_slice());
        })
    });
}

criterion_group!(
    name = benches;
    //config = Criterion::default().sample_size(10);
    config = Criterion::default();
    targets = lagrange_group);
criterion_main!(benches);