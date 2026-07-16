from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

from type_defs import (
    ContextT,
    GmwfsContextT,
    GmwfsPublicKeyT,
    GmwfsSecretKeyT,
    GmwfsSignatureT,
    GroupElementT,
    PublicKeyT,
    SecretKeyT,
    SetElementT,
    SignatureT,
)


class SignatureScheme(Protocol[PublicKeyT, SecretKeyT, SignatureT, ContextT]):
    name: str
    context: ContextT

    # Core API
    def keygen(self) -> tuple[PublicKeyT, SecretKeyT, ContextT]: ...

    def sign(
        self,
        secret_key: SecretKeyT,
        public_key: PublicKeyT,
        message: bytes,
        context: ContextT,
    ) -> tuple[SignatureT, ContextT]: ...

    def verify(
        self,
        public_key: PublicKeyT,
        message: bytes,
        signature: SignatureT,
        context: ContextT,
    ) -> bool: ...

    # Serialized, end-to-end API
    def keygen_bytes(self) -> tuple[bytes, bytes]: ...

    def sign_bytes(
        self,
        secret_key: bytes,
        public_key: bytes,
        message: bytes,
    ) -> bytes: ...

    def verify_bytes(
        self,
        public_key: bytes,
        message: bytes,
        signature: bytes,
    ) -> bool: ...


type SignatureSchemeT = SignatureScheme

type GmwfsScheme = SignatureScheme


@dataclass
class GmwfsAdapter[SetElementT]:
    name: str
    scheme: GmwfsScheme
    context: GmwfsContextT = None

    def keygen(
        self,
    ) -> tuple[
        GmwfsPublicKeyT[SetElementT],
        GmwfsSecretKeyT[GroupElementT],
        GmwfsContextT,
    ]:
        return self.scheme.keygen()

    def sign(
        self,
        secret_key: GmwfsSecretKeyT[GroupElementT],
        public_key: GmwfsPublicKeyT[SetElementT],
        message: bytes,
        context: GmwfsContextT,
    ) -> tuple[GmwfsSignatureT[GroupElementT], GmwfsContextT]:
        signature = self.scheme.sign(
            secret_key,
            public_key,
            message,
            context,
        )

        return signature, context

    def verify(
        self,
        public_key: GmwfsPublicKeyT[SetElementT],
        message: bytes,
        signature: GmwfsSignatureT[GroupElementT],
        context: GmwfsContextT,
    ) -> bool:
        return bool(
            self.scheme.verify(
                public_key,
                message,
                signature,
                context,
            )
        )

    def keygen_bytes(self) -> tuple[bytes, bytes]:
        return self.scheme.keygen_bytes()

    def sign_bytes(
        self,
        secret_key: bytes,
        public_key: bytes,
        message: bytes,
    ) -> bytes:
        return self.scheme.sign_bytes(
            secret_key,
            public_key,
            message,
        )

    def verify_bytes(
        self,
        public_key: bytes,
        message: bytes,
        signature: bytes,
    ) -> bool:
        return bool(
            self.scheme.verify_bytes(
                public_key,
                message,
                signature,
            )
        )


def load_scheme(name: str) -> SignatureSchemeT:
    normalized = name.lower().replace("_", "-")

    if normalized in {"owl-mlip", "mlip"}:
        from owl_mlip.main import owl_mlip

        return GmwfsAdapter(name="OWL-mLIP", scheme=owl_mlip)

    if normalized in {"owl-lip", "lip"}:
        from owl_lip.main import owl_lip

        return GmwfsAdapter(name="OWL-LIP", scheme=owl_lip)

    raise ValueError(f"Unknown scheme {name!r}. Choose one of: owl-mlip, owl-lip.")
