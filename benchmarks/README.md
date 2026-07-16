# OWL-GMW-FS Benchmark Suite

This directory contains reproducible benchmarks for OWL-mLIP and OWL-LIP.

## Contents

```text
benchmarks/
├── benchmark.py
├── benchmark_keygen.py
├── benchmark_sign.py
├── benchmark_verify.py
├── benchmark_sizes.py
├── benchmark_scaling.py
├── profile.py
├── common.py
├── scheme_adapter.py
└── results/
    └── plots/
```

## Important setup

The benchmark code assumes the following OWL-mLIP API:

```python
pk, sk, orbit = owl_mLIP.keygen()
signature = owl_mLIP.sign(sk, orbit, message)
valid = owl_mLIP.verify(pk, orbit, message, signature)
```

Scheme-specific imports and serialization are isolated in `scheme_adapter.py`.

Before publishing size results, replace the `pickle` fallback with the canonical
public-key, secret-key, and signature encoders used by the protocol. Pickle size
is an implementation-storage measurement, not a cryptographic wire-format size.

## Run all benchmarks

From the repository root:

```bash
python benchmarks/benchmark.py \
    --scheme owl-mlip \
    --iterations 30 \
    --warmup 3 \
    --message-size 32 \
    --output-dir benchmarks/results/owl_mlip
```

## Run individual benchmarks

```bash
python benchmarks/benchmark_keygen.py --scheme owl-mlip
python benchmarks/benchmark_sign.py --scheme owl-mlip
python benchmarks/benchmark_verify.py --scheme owl-mlip
python benchmarks/benchmark_sizes.py --scheme owl-mlip
```

Use `--scheme owl-lip` for OWL-LIP.

## Scaling benchmark

Example for `LOGN`:

```bash
python benchmarks/benchmark_scaling.py \
    --scheme owl-mlip \
    --parameter-module owl_mLIP.params \
    --parameter LOGN \
    --values 8 9 10 \
    --iterations 10
```

Parameter mutation only works when the implementation reads parameters
dynamically. When constants are imported directly into other modules, run each
parameter configuration in a fresh Python process.

## Profiling

```bash
python benchmarks/profile.py \
    --scheme owl-mlip \
    --iterations 3 \
    --output benchmarks/results/owl_mlip/profile.prof
```

Inspect the profile interactively with tools such as `snakeviz`:

```bash
python -m pip install snakeviz
snakeviz benchmarks/results/owl_mlip/profile.prof
```

## Output files

A complete run produces:

```text
raw.csv
sizes.csv
summary.json
```

`raw.csv` preserves every timing sample. `summary.json` contains descriptive
statistics and environment metadata.

## Methodology

- Timings use `time.perf_counter_ns()`.
- Warm-up runs are excluded.
- Signing uses a fresh random message for every iteration.
- Every generated signature is verified.
- Key generation, signing, and verification are measured separately.
- Raw measurements are retained instead of reporting only averages.

For publication-quality results, also record CPU model, operating-system
version, Python version, dependency versions, Git commit, parameter set, and
whether CPU frequency scaling or background processes were controlled.

## Recommended reporting

Report at least:

| Scheme | PK bytes | SK bytes | Signature bytes | KeyGen mean | Sign mean | Verify mean |
|---|---:|---:|---:|---:|---:|---:|

Also report the median, standard deviation, number of repetitions, machine
configuration, exact Git commit, and protocol parameter set.

## Security note

This benchmark validates functional correctness only. It does not evaluate
side-channel resistance, constant-time behavior, entropy quality, or the
cryptographic security of the construction.
