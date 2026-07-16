from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Any

from benchmarks.benchmark_keygen import benchmark_keygen
from benchmarks.benchmark_sign import benchmark_sign
from benchmarks.benchmark_verify import benchmark_verify
from benchmarks.common import environment_info, write_json, write_raw_csv
from benchmarks.scheme_adapter import load_scheme


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the complete signature benchmark."
    )
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
    args.output_dir.mkdir(parents=True, exist_ok=True)

    keygen_rows, keygen_summary = benchmark_keygen(
        scheme,
        args.iterations,
        args.warmup,
    )
    sign_rows, sign_summary = benchmark_sign(
        scheme,
        args.iterations,
        args.warmup,
        args.message_size,
    )
    verify_rows, verify_summary = benchmark_verify(
        scheme,
        args.iterations,
        args.warmup,
        args.message_size,
    )

    # Final correctness check using the core API.
    public_key, secret_key, context = scheme.keygen()
    message = os.urandom(args.message_size)

    signature, signing_context = scheme.sign(
        secret_key,
        public_key,
        message,
        context,
    )

    if not scheme.verify(
        public_key,
        message,
        signature,
        signing_context,
    ):
        raise RuntimeError("Final core correctness check failed.")

    # Measure serialized sizes using the byte API.
    public_key_bytes, secret_key_bytes = scheme.keygen_bytes()

    signature_bytes = scheme.sign_bytes(
        secret_key_bytes,
        public_key_bytes,
        message,
    )

    if not scheme.verify_bytes(
        public_key_bytes,
        message,
        signature_bytes,
    ):
        raise RuntimeError("Final byte API correctness check failed.")

    size_row: dict[str, str | int] = {
        "scheme": scheme.name,
        "public_key_bytes": len(public_key_bytes),
        "secret_key_bytes": len(secret_key_bytes),
        "signature_bytes": len(signature_bytes),
        "message_bytes": len(message),
    }

    write_raw_csv(
        args.output_dir / "raw.csv",
        [*keygen_rows, *sign_rows, *verify_rows],
    )
    write_raw_csv(
        args.output_dir / "sizes.csv",
        [size_row],
    )

    summary: dict[str, Any] = {
        "environment": environment_info(),
        "configuration": {
            "scheme": scheme.name,
            "iterations": args.iterations,
            "warmup": args.warmup,
            "message_size_bytes": args.message_size,
        },
        "timings": {
            "keygen_core": keygen_summary,
            "sign_core": sign_summary,
            "verify_core": verify_summary,
        },
        "sizes": size_row,
    }

    write_json(
        args.output_dir / "summary.json",
        summary,
    )

    print(f"Benchmark complete. Results written to {args.output_dir}/")


if __name__ == "__main__":
    main()
