# GMW-FS
# A post-quantum digital signature protocol using cryptographic group actions

from collections.abc import Callable
from dataclasses import dataclass, field
import os

from cga import Cga
from gmwfs_codec import GmwfsEncoder
from rngcontext import RngContext, make_rng
from type_defs import (
    GmwfsContextT,
    GmwfsPublicKeyT,
    GmwfsSecretKeyT,
    GmwfsSignatureT,
    GroupElementT,
    SetElementT,
)


@dataclass
class Gmwfs[GroupElementT, SetElementT]:
    name: str
    cga: Cga[GroupElementT, SetElementT]
    c: int
    k: int
    r: int
    hash: Callable[[bytes], bytes]
    expand: Callable[[bytes], tuple[list[int], list[int], list[int]]]
    chg_size: int
    context: SetElementT | None
    rngcontext: RngContext = field(init=False)

    def __post_init__(self):
        self.encoder = GmwfsEncoder(self.cga, self.chg_size)
        x = int.from_bytes(os.urandom(1), byteorder="big")
        self.rngcontext = make_rng(x)

    def keygen(
        self,
    ) -> tuple[
        GmwfsPublicKeyT[SetElementT],
        GmwfsSecretKeyT[GroupElementT],
        GmwfsContextT,
    ]:
        orbit = (
            self.context
            if self.context is not None
            else self.cga.set.sample(self.rngcontext)
        )
        secret_key = [self.cga.group.sample(self.rngcontext) for i in range(self.c)]
        public_key = [
            self.cga.action(orbit, self.cga.group.inverse(secret_key[i]))
            for i in range(self.c)
        ]

        return (orbit, public_key), secret_key, None

    def sign(
        self,
        secret_key: GmwfsSecretKeyT[GroupElementT],
        public_key: GmwfsPublicKeyT[SetElementT],
        message: bytes,
        context: GmwfsContextT = None,
    ) -> tuple[GmwfsSignatureT[GroupElementT], GmwfsContextT]:
        orbit, _ = public_key
        h_i = [self.cga.group.sample(self.rngcontext) for i in range(self.r)]
        t_i = [self.cga.action(orbit, h_i[i]) for i in range(self.r)]

        hash_input = message
        for i in t_i:
            hash_input += self.cga.set.encode(i)

        chg = self.hash(hash_input)
        chg_c, chg_nc, chg_val = self.expand(chg)

        f_i = [self.cga.group.placeholder()] * self.r
        for k in range(self.k):
            f_i[chg_nc[k]] = self.cga.group.operate(
                h_i[chg_nc[k]], secret_key[chg_val[k]]
            )

        for k in range(self.r - self.k):
            f_i[chg_c[k]] = h_i[chg_c[k]]

        return ((chg, f_i), context)

    def verify(
        self,
        public_key: GmwfsPublicKeyT[SetElementT],
        message: bytes,
        signature: GmwfsSignatureT[GroupElementT],
        context: GmwfsContextT = None,
    ) -> bool:
        orbit, public = public_key
        chg, f_i = signature
        chg_c, chg_nc, chg_val = self.expand(chg)
        t_i = [self.cga.set.placeholder()] * self.r
        for k in range(self.k):
            t_i[chg_nc[k]] = self.cga.action(public[chg_val[k]], f_i[chg_nc[k]])

        for k in range(self.r - self.k):
            t_i[chg_c[k]] = self.cga.action(orbit, f_i[chg_c[k]])

        hash_input = message
        for i in t_i:
            hash_input += self.cga.set.encode(i)
        chg_2 = self.hash(hash_input)

        return chg == chg_2

    def keygen_bytes(self) -> tuple[bytes, bytes]:
        print("Generating Key...")
        public_key, secret_key, _ = self.keygen()

        pk = self.encoder.encode_public_key(public_key)
        sk = self.encoder.encode_secret_key(secret_key)

        print("Key Generated!")
        return pk, sk

    def sign_bytes(self, sk: bytes, pk: bytes, message: bytes) -> bytes:
        print("Signing Message...")
        secret_key = self.encoder.decode_secret_key(sk)
        public_key = self.encoder.decode_public_key(pk)

        sign, _ = self.sign(secret_key, public_key, message)

        print("Message Signed!")
        return self.encoder.encode_signature(sign)

    def verify_bytes(self, pk: bytes, message: bytes, sign: bytes) -> bool:
        print("Verifying Message...")
        public_key = self.encoder.decode_public_key(pk)
        decoded_sign = self.encoder.decode_signature(sign)

        valid_sign = self.verify(public_key, message, decoded_sign)

        if valid_sign:
            print("Verification Success")
            return True
        print("Verification Failed")
        return False

    def testing(self):
        pk, sk = self.keygen_bytes()

        msg = os.urandom(32)

        sig = self.sign_bytes(sk, pk, msg)

        self.verify_bytes(pk, msg, sig)
