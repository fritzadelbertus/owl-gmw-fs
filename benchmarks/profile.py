from __future__ import annotations

import argparse
import cProfile
import os
from pathlib import Path
import pstats

from benchmarks.scheme_adapter import load_scheme


def workload(scheme_name: str, iterations: int, message_size: int) -> None:
    scheme = load_scheme(scheme_name)

    for _ in range(iterations):
        public_key, secret_key, orbit = scheme.keygen()
        message = os.urandom(message_size)
        signature = scheme.sign(secret_key, public_key, message, orbit)
        if not scheme.verify(public_key, orbit, message, signature):
            raise RuntimeError("Verification failed during profiling.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Profile the complete signature workflow."
    )
    parser.add_argument("--scheme", default="owl-mlip")
    parser.add_argument("--iterations", type=int, default=3)
    parser.add_argument("--message-size", type=int, default=32)
    parser.add_argument("--output", type=Path, default=Path("results/profile.prof"))
    parser.add_argument("--sort", default="cumulative")
    parser.add_argument("--top", type=int, default=40)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)

    profiler = cProfile.Profile()
    profiler.enable()
    workload(args.scheme, args.iterations, args.message_size)
    profiler.disable()
    profiler.dump_stats(args.output)

    stats = pstats.Stats(profiler).strip_dirs().sort_stats(args.sort)
    stats.print_stats(args.top)
    print(f"\nProfile saved to {args.output}")


if __name__ == "__main__":
    main()
