from cga import Cga
from type_defs import (
    GmwfsPublicKeyT,
    GmwfsSecretKeyT,
    GmwfsSignatureT,
    GroupElementT,
    SetElementT,
)


class GmwfsEncoder[GroupElementT, SetElementT]:
    def __init__(self, cga: Cga[GroupElementT, SetElementT], chg_size: int):
        self.cga = cga
        self.chg_size = chg_size

    def encode_signature(self, sign: GmwfsSignatureT[GroupElementT]) -> bytes:
        chg, f_i = sign
        encoded_sign = chg
        for f in f_i:
            encoded_sign += self.cga.group.encode(f)
        return encoded_sign

    def decode_signature(self, sign: bytes) -> GmwfsSignatureT[GroupElementT]:
        chg = sign[: self.chg_size]
        right_part = sign[self.chg_size :]
        f_i_byte = [
            right_part[i : i + self.cga.group.element_length()]
            for i in range(0, len(right_part), self.cga.group.element_length())
        ]

        f_i = [self.cga.group.decode(f) for f in f_i_byte]

        return chg, f_i

    def encode_public_key(self, public_key: GmwfsPublicKeyT[SetElementT]) -> bytes:
        orbit, public = public_key
        pk = self.cga.set.encode(orbit)
        for key in public:
            pk += self.cga.set.encode(key)
        return pk

    def decode_public_key(self, pk: bytes) -> GmwfsPublicKeyT[SetElementT]:
        length = self.cga.set.element_length()
        orb = pk[:length]
        public_key_part = pk[length:]
        public_key_bytes = [
            public_key_part[j : j + length]
            for j in range(0, len(public_key_part), length)
        ]
        public_key = [self.cga.set.decode(i) for i in public_key_bytes]
        orbit = self.cga.set.decode(orb)
        return orbit, public_key

    def encode_secret_key(self, secret_key: GmwfsSecretKeyT[GroupElementT]) -> bytes:
        sk = b""
        for key in secret_key:
            sk += self.cga.group.encode(key)
        return sk

    def decode_secret_key(self, sk: bytes) -> GmwfsSecretKeyT[GroupElementT]:
        length = self.cga.group.element_length()
        return [
            self.cga.group.decode(i)
            for i in [sk[j : j + length] for j in range(0, len(sk), length)]
        ]
