from __future__ import annotations

import argparse
from functools import partial
import os
from pathlib import Path
from typing import Any

from benchmarks.common import (
    environment_info,
    measure_ns,
    print_summary,
    summarize_ns,
    write_json,
    write_raw_csv,
)
from benchmarks.scheme_adapter import SignatureSchemeT, load_scheme


def benchmark_verify(
    scheme: SignatureSchemeT,
    iterations: int,
    warmup: int,
    message_size: int,
) -> tuple[list[dict[str, Any]], dict[str, float | int]]:
    public_key, secret_key, context = scheme.keygen()

    samples: list[tuple[bytes, Any, Any]] = []

    for _ in range(iterations + warmup):
        message = os.urandom(message_size)

        signature, signing_context = scheme.sign(
            secret_key,
            public_key,
            message,
            context,
        )

        samples.append(
            (
                message,
                signature,
                signing_context,
            )
        )

    for message, signature, signing_context in samples[:warmup]:
        verified = scheme.verify(
            public_key,
            message,
            signature,
            signing_context,
        )

        if not verified:
            raise RuntimeError("Verification failed during warm-up.")

    rows: list[dict[str, Any]] = []
    timings_ns: list[int] = []

    for iteration, sample in enumerate(samples[warmup:]):
        message, signature, signing_context = sample

        verify_call = partial(
            scheme.verify,
            public_key,
            message,
            signature,
            signing_context,
        )

        verified, elapsed_ns = measure_ns(verify_call)

        if not verified:
            raise RuntimeError(f"Verification failed at iteration {iteration}.")

        timings_ns.append(elapsed_ns)
        rows.append(
            {
                "scheme": scheme.name,
                "operation": "verify_core",
                "iteration": iteration,
                "message_size_bytes": message_size,
                "elapsed_ns": elapsed_ns,
                "elapsed_ms": elapsed_ns / 1_000_000,
            }
        )

    return rows, summarize_ns(timings_ns)


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark signature verification.")
    parser.add_argument("--scheme", default="owl-mlip")
    parser.add_argument("--iterations", type=int, default=30)
    parser.add_argument("--warmup", type=int, default=3)
    parser.add_argument("--message-size", type=int, default=32)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results"),
    )
    args = parser.parse_args()

    scheme = load_scheme(args.scheme)

    rows, summary = benchmark_verify(
        scheme,
        args.iterations,
        args.warmup,
        args.message_size,
    )

    write_raw_csv(
        args.output_dir / "verify_raw.csv",
        rows,
    )
    write_json(
        args.output_dir / "verify_summary.json",
        {
            "environment": environment_info(),
            "scheme": scheme.name,
            "operation": "verify_core",
            "message_size_bytes": args.message_size,
            "summary": summary,
        },
    )
    print_summary("Core verification", summary)


if __name__ == "__main__":
    main()
