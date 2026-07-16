import hashlib

import numpy as np
from numpy.typing import NDArray


class RngContext:
    """
    This class is used to generate deterministically randomness.

    This is equivalent to the rngcontext of the C reference implementation.
    """

    def __init__(self, seed: NDArray[np.uint8]) -> None:
        self.shake256 = hashlib.shake_256()
        self.shake256.update(seed.tobytes())
        self.i = 0

    def random(self, size: int) -> NDArray[np.uint8]:
        r = np.frombuffer(self.shake256.digest(self.i + size), dtype=np.uint8)
        r = r[self.i : self.i + size]
        self.i += size
        return r


def make_rng(seed: int) -> RngContext:
    return RngContext(np.array([seed], dtype=np.uint8))


def shake256x4(message: bytes, num: int) -> list[int]:
    """
    SHAKE256x4 (see Section 3.3.2)

    Inputs:
        - message: input to SHAKE256
        - num: number of 64-bit output elements requested

    Outputs:
        - y : array of num 64-bit elements
    """
    digest: list[bytes] = [b""] * 4
    for i in range(4):
        shake = hashlib.shake_256()
        shake.update(message + bytes([i]))
        digest[i] = shake.digest(num * 8 // 4)
        # *8 (8 bytes in 64 bits) / 4 (for the 4 streams)

    # Convert the four streams in the interleaved form
    j = 0
    y: list[int] = [0] * int(num)
    for i in range(int(num / 4)):
        y[j + 0] = int.from_bytes(digest[0][i * 8 : (i + 1) * 8], byteorder="little")
        y[j + 1] = int.from_bytes(digest[1][i * 8 : (i + 1) * 8], byteorder="little")
        y[j + 2] = int.from_bytes(digest[2][i * 8 : (i + 1) * 8], byteorder="little")
        y[j + 3] = int.from_bytes(digest[3][i * 8 : (i + 1) * 8], byteorder="little")
        j = j + 4
    return y
