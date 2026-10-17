"""Read the exact AO/MO data shared by the PySCF and C++ CASSCF runners."""

import hashlib
import json
from pathlib import Path
import struct
import sys

import numpy as np

MAGIC = b"XMCASAO1"
HEADER = struct.Struct("<8sIIIIId")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read(fixture, metadata):
    fixture, metadata = Path(fixture), Path(metadata)
    data = json.loads(metadata.read_text())
    if sha256(fixture) != data["fixture_sha256"]:
        raise ValueError("fixture SHA-256 does not match metadata")
    with fixture.open("rb") as stream:
        magic, nbf, ncore, ncas, nelecas, spin, nuclear = HEADER.unpack(
            stream.read(HEADER.size))
    if magic != MAGIC or nbf < 1 or ncore + ncas > nbf:
        raise ValueError("invalid XMCASAO1 header")
    npair = nbf * (nbf + 1) // 2
    neri = npair * (npair + 1) // 2
    expected = HEADER.size + 3 * nbf * nbf * 8 + neri * 8
    if fixture.stat().st_size != expected:
        raise ValueError("invalid XMCASAO1 byte count")
    if sys.byteorder != "little":
        raise ValueError("XMCASAO1 requires a little-endian host")
    shape = (nbf, nbf)
    matrix_bytes = nbf * nbf * 8
    s = np.memmap(fixture, dtype="<f8", mode="r", offset=HEADER.size,
                  shape=shape, order="C")
    h = np.memmap(fixture, dtype="<f8", mode="r", offset=HEADER.size + matrix_bytes,
                  shape=shape, order="C")
    u = np.memmap(fixture, dtype="<f8", mode="r", offset=HEADER.size + 2 * matrix_bytes,
                  shape=shape, order="C")
    eri = np.memmap(fixture, dtype="<f8", mode="r",
                    offset=HEADER.size + 3 * matrix_bytes, shape=(neri,))
    if max(np.max(np.abs(s - s.T)), np.max(np.abs(h - h.T)),
           np.max(np.abs(u.T @ s @ u - np.eye(nbf)))) > 1e-8:
        raise ValueError("invalid AO symmetry or initial MO orthonormality")
    return dict(metadata=data, nbf=nbf, ncore=ncore, ncas=ncas,
                nelecas=nelecas, spin=spin, nuclear=nuclear,
                overlap=s, hcore=h, initial_mo=u, eri=eri)
