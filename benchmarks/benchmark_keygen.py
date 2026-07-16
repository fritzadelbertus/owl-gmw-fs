from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from benchmarks.common import (
    environment_info,
    measure_ns,
    print_summary,
    summarize_ns,
    warm_up,
    write_json,
    write_raw_csv,
)
from benchmarks.scheme_adapter import SignatureSchemeT, load_scheme


def benchmark_keygen(
    scheme: SignatureSchemeT,
    iterations: int,
    warmup: int,
) -> tuple[list[dict[str, Any]], dict[str, float | int]]:
    warm_up(scheme.keygen, warmup)

    rows: list[dict[str, Any]] = []
    timings_ns: list[int] = []

    for iteration in range(iterations):
        _, elapsed_ns = measure_ns(scheme.keygen)

        timings_ns.append(elapsed_ns)
        rows.append(
            {
                "scheme": scheme.name,
                "operation": "keygen",
                "iteration": iteration,
                "elapsed_ns": elapsed_ns,
                "elapsed_ms": elapsed_ns / 1_000_000,
            }
        )

    return rows, summarize_ns(timings_ns)


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark signature key generation.")
    parser.add_argument("--scheme", default="owl-mlip")
    parser.add_argument("--iterations", type=int, default=30)
    parser.add_argument("--warmup", type=int, default=3)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results"),
    )
    args = parser.parse_args()

    scheme = load_scheme(args.scheme)

    rows, summary = benchmark_keygen(
        scheme,
        args.iterations,
        args.warmup,
    )

    write_raw_csv(
        args.output_dir / "keygen_raw.csv",
        rows,
    )
    write_json(
        args.output_dir / "keygen_summary.json",
        {
            "environment": environment_info(),
            "scheme": scheme.name,
            "summary": summary,
        },
    )
    print_summary("Key generation", summary)


if __name__ == "__main__":
    main()
