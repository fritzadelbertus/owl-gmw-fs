from __future__ import annotations

from collections.abc import Callable, Iterable
import csv
from dataclasses import dataclass
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time
from typing import Any


@dataclass(frozen=True)
class TimingResult:
    operation: str
    iteration: int
    elapsed_ns: int

    @property
    def elapsed_ms(self) -> float:
        return self.elapsed_ns / 1_000_000


def measure_ns(func: Callable[[], Any]) -> tuple[Any, int]:
    start = time.perf_counter_ns()
    result = func()
    elapsed_ns = time.perf_counter_ns() - start
    return result, elapsed_ns


def warm_up(func: Callable[[], Any], iterations: int) -> None:
    for _ in range(iterations):
        func()


def summarize_ns(values_ns: Iterable[int]) -> dict[str, float | int]:
    values = list(values_ns)
    if not values:
        raise ValueError("Cannot summarize an empty timing sequence.")

    values_ms = [value / 1_000_000 for value in values]
    mean_ms = statistics.fmean(values_ms)

    return {
        "runs": len(values_ms),
        "mean_ms": mean_ms,
        "median_ms": statistics.median(values_ms),
        "stdev_ms": statistics.stdev(values_ms) if len(values_ms) > 1 else 0.0,
        "min_ms": min(values_ms),
        "max_ms": max(values_ms),
        "operations_per_second": 1000.0 / mean_ms if mean_ms > 0 else float("inf"),
    }


def write_raw_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("No rows were provided.")

    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)

    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as file:
        json.dump(data, file, indent=2, sort_keys=True)


def environment_info() -> dict[str, str]:
    return {
        "platform": platform.platform(),
        "python": sys.version.replace("\n", " "),
        "processor": platform.processor() or "unknown",
        "machine": platform.machine(),
        "cpu_count": str(os.cpu_count()),
        "git_commit": current_git_commit(),
    }


def current_git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except OSError, subprocess.CalledProcessError:
        return "unknown"


def print_summary(operation: str, summary: dict[str, float | int]) -> None:
    print(f"\n{operation}")
    print("-" * len(operation))
    print(f"runs:       {summary['runs']}")
    print(f"mean:       {summary['mean_ms']:.6f} ms")
    print(f"median:     {summary['median_ms']:.6f} ms")
    print(f"stdev:      {summary['stdev_ms']:.6f} ms")
    print(f"min:        {summary['min_ms']:.6f} ms")
    print(f"max:        {summary['max_ms']:.6f} ms")
    print(f"throughput: {summary['operations_per_second']:.3f} ops/s")
