# AGENTS.md

# OWL-GMW-FS Contributor Guide

This project prioritizes **correctness, readability, reproducibility, and strong typing**
over writing the shortest possible code.

The codebase should be understandable without requiring readers to inspect every
implementation detail.

---

# Core Principles

## 1. Everything must be typed

Every public function, method, class attribute, and dataclass field must have
explicit type annotations.

Do not leave return types implicit.

Prefer

```python
def group_inverse(g: GroupElementT) -> GroupElementT:
```

instead of

```python
def group_inverse(g):
```

---

## 2. The function signature should explain the function

A reader should understand what a function does by reading only

- its name
- its parameters
- its return type

Avoid comments that explain *what* the code does.

Instead, improve

- naming
- typing
- decomposition into smaller functions.

Example:

Good

```python
def canonical_pm(matrix: GroupElementT) -> GroupElementT:
```

Bad

```python
def f(x):
    # Normalize the matrix into its canonical representative.
```

---

## 3. Comments explain **why**, not **what**

Comments should only exist when they explain reasoning that cannot easily be
expressed in code.

Avoid comments such as

```python
# Increment i
i += 1
```

Acceptable comments include

```python
# The determinant is guaranteed to be ±1 because every generator is unimodular.
```

Function comments, if necessary, should be at most **one line immediately below
the function definition.**

---

## 4. Small functions

Functions should perform one conceptual task.

If a function requires several paragraphs of comments, split it into multiple
functions.

---

## 5. Prefer explicit code over clever code

Readable code is preferred over compact code.

Avoid nested comprehensions or deeply nested expressions when they reduce
clarity.

---

## 6. Deterministic behaviour

Randomness must come from `RngContext`.

Do not use

- `random`
- `numpy.random`
- `os.urandom`

unless the function explicitly generates an external seed.

---

## 7. Cryptographic correctness first

Never introduce optimizations that change mathematical behaviour.

Correct finite-field arithmetic is more important than speed.

Avoid relying on NumPy overflow semantics unless the algorithm explicitly
requires arithmetic modulo `2^n`.

---

## 8. Keep abstractions consistent

Each abstraction should expose a single mathematical concept.

For example

- `Group`
- `Set`
- `GroupAction`

should not leak implementation-specific details.

---

## 9. Prefer descriptive names

Names should describe mathematical meaning.

Prefer

```python
group_inverse
group_action
sample_group
```

over

```python
ginv
act
sample
```

Temporary variables may be short when they follow mathematical notation.

---

## 10. Avoid unnecessary comments

If removing a comment makes the code harder to understand,
rename variables or split the function instead.

---

# Style

Formatting is enforced by Ruff.

Imports are sorted automatically.

Typing is checked using mypy.

Do not disable lint rules unless there is a documented mathematical reason.

---

# Testing

Every public function should have tests whenever practical.

Important properties include

- inverse correctness
- encode/decode round trips
- group laws
- deterministic RNG behaviour
- serialization compatibility

Property-based tests are preferred where applicable.

---

# Pull Requests

Changes should

- preserve deterministic behaviour
- preserve typing
- pass Ruff
- pass mypy
- pass pytest
- maintain mathematical correctness
