from __future__ import annotations

import argparse
from functools import partial
import importlib
import os
from pathlib import Path
from typing import Any

from benchmarks.common import measure_ns, write_raw_csv
from benchmarks.scheme_adapter import load_scheme


def set_parameter(
    module_name: str,
    parameter: str,
    value: int,
) -> None:
    """
    Mutate a parameter module before executing a benchmark point.

    This only works when the scheme reads parameters dynamically. If constants
    are imported using `from params import X`, execute each parameter point in
    a fresh subprocess instead.
    """
    module = importlib.import_module(module_name)

    if not hasattr(module, parameter):
        raise AttributeError(f"{module_name} has no parameter {parameter!r}.")

    setattr(module, parameter, value)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Benchmark scaling over one parameter."
    )
    parser.add_argument("--scheme", default="owl-mlip")
    parser.add_argument("--parameter-module", required=True)
    parser.add_argument("--parameter", required=True)
    parser.add_argument(
        "--values",
        type=int,
        nargs="+",
        required=True,
    )
    parser.add_argument("--iterations", type=int, default=10)
    parser.add_argument("--message-size", type=int, default=32)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results"),
    )
    args = parser.parse_args()

    rows: list[dict[str, Any]] = []

    for value in args.values:
        set_parameter(
            args.parameter_module,
            args.parameter,
            value,
        )

        scheme = load_scheme(args.scheme)

        for iteration in range(args.iterations):
            keygen_result, keygen_ns = measure_ns(scheme.keygen)
            public_key, secret_key, context = keygen_result

            message = os.urandom(args.message_size)

            sign_call = partial(
                scheme.sign,
                secret_key,
                public_key,
                message,
                context,
            )

            sign_result, sign_ns = measure_ns(sign_call)
            signature, signing_context = sign_result

            verify_call = partial(
                scheme.verify,
                public_key,
                message,
                signature,
                signing_context,
            )

            verified, verify_ns = measure_ns(verify_call)

            if not verified:
                raise RuntimeError(
                    f"Verification failed for "
                    f"{args.parameter}={value}, "
                    f"iteration={iteration}."
                )

            rows.append(
                {
                    "scheme": scheme.name,
                    "parameter": args.parameter,
                    "parameter_value": value,
                    "iteration": iteration,
                    "message_size_bytes": args.message_size,
                    "keygen_ns": keygen_ns,
                    "sign_ns": sign_ns,
                    "verify_ns": verify_ns,
                    "keygen_ms": keygen_ns / 1_000_000,
                    "sign_ms": sign_ns / 1_000_000,
                    "verify_ms": verify_ns / 1_000_000,
                }
            )

    output = args.output_dir / f"scaling_{args.parameter.lower()}.csv"

    write_raw_csv(output, rows)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
