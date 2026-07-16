from __future__ import annotations

import argparse
import os
from pathlib import Path

from benchmarks.common import (
    environment_info,
    write_json,
    write_raw_csv,
)
from benchmarks.scheme_adapter import load_scheme


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Measure serialized key and signature sizes."
    )
    parser.add_argument("--scheme", default="owl-mlip")
    parser.add_argument("--message-size", type=int, default=32)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results"),
    )
    args = parser.parse_args()

    scheme = load_scheme(args.scheme)

    public_key_bytes, secret_key_bytes = scheme.keygen_bytes()

    message = os.urandom(args.message_size)

    signature_bytes = scheme.sign_bytes(
        secret_key_bytes,
        public_key_bytes,
        message,
    )

    verified = scheme.verify_bytes(
        public_key_bytes,
        message,
        signature_bytes,
    )

    if not verified:
        raise RuntimeError("The generated serialized signature did not verify.")

    sizes: dict[str, str | int] = {
        "scheme": scheme.name,
        "public_key_bytes": len(public_key_bytes),
        "secret_key_bytes": len(secret_key_bytes),
        "signature_bytes": len(signature_bytes),
        "message_bytes": len(message),
    }

    write_raw_csv(
        args.output_dir / "sizes.csv",
        [sizes],
    )
    write_json(
        args.output_dir / "sizes.json",
        {
            "environment": environment_info(),
            "sizes": sizes,
        },
    )

    for key, value in sizes.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
