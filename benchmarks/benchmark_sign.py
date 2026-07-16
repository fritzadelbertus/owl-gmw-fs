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
    warm_up,
    write_json,
    write_raw_csv,
)
from benchmarks.scheme_adapter import SignatureSchemeT, load_scheme


def benchmark_sign(
    scheme: SignatureSchemeT,
    iterations: int,
    warmup: int,
    message_size: int,
) -> tuple[list[dict[str, Any]], dict[str, float | int]]:
    public_key, secret_key, context = scheme.keygen()

    def one_sign() -> tuple[Any, Any]:
        message = os.urandom(message_size)

        return scheme.sign(
            secret_key,
            public_key,
            message,
            context,
        )

    warm_up(one_sign, warmup)

    rows: list[dict[str, Any]] = []
    timings_ns: list[int] = []

    for iteration in range(iterations):
        message = os.urandom(message_size)

        sign_call = partial(
            scheme.sign,
            secret_key,
            public_key,
            message,
            context,
        )

        result, elapsed_ns = measure_ns(sign_call)
        signature, signing_context = result

        valid = scheme.verify(
            public_key,
            message,
            signature,
            signing_context,
        )

        if not valid:
            raise RuntimeError(f"Verification failed at iteration {iteration}.")

        timings_ns.append(elapsed_ns)
        rows.append(
            {
                "scheme": scheme.name,
                "operation": "sign_core",
                "iteration": iteration,
                "message_size_bytes": message_size,
                "elapsed_ns": elapsed_ns,
                "elapsed_ms": elapsed_ns / 1_000_000,
            }
        )

    return rows, summarize_ns(timings_ns)


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark signature generation.")
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

    rows, summary = benchmark_sign(
        scheme,
        args.iterations,
        args.warmup,
        args.message_size,
    )

    write_raw_csv(
        args.output_dir / "sign_raw.csv",
        rows,
    )
    write_json(
        args.output_dir / "sign_summary.json",
        {
            "environment": environment_info(),
            "scheme": scheme.name,
            "operation": "sign_core",
            "message_size_bytes": args.message_size,
            "summary": summary,
        },
    )
    print_summary("Core signing", summary)


if __name__ == "__main__":
    main()
