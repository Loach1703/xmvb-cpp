#!/usr/bin/env python3
"""Independent C++ TNVB protocol: singlet SU(2), fixed D, DIRECT, FULL only.

The paired core is an owned copy of stable worker revision
fbe78374bdcf1dd73930fe5f5c3f0d403d5c4cbf6e30d093b9a8a16fa6cc3e73.
Energy requests never invoke the total-state/RDM branch of that core.
"""

import os
import struct
import sys
import time
import traceback

# This process (and its spawned workers) owns the numerical runtime.
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["GOTO_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["XMVB_BLOCK2_RDM_CHECK_TOL"] = "1e-7"
os.environ.pop("XMVB_BLOCK2_SUM_ENERGY_TOL", None)
os.environ.pop("XMVB_BLOCK2_VALIDATE_DIRECT_D3G", None)

import numpy as np
import scipy.linalg
import pair_core as core

REQUEST_MAGIC = b"CPTNREQ1"
RESPONSE_MAGIC = b"CPTNRES1"
VERSION = 1


def read_array(stream, count, dtype, shape=None):
    itemsize = np.dtype(dtype).itemsize
    value = np.frombuffer(core.read_exact(stream, count * itemsize), dtype=dtype)
    return value.copy() if shape is None else value.reshape(shape, order="F").copy()


def read_request(path):
    with open(path, "rb") as stream:
        if core.read_exact(stream, 8) != REQUEST_MAGIC:
            raise ValueError("invalid C++ TNVB request magic")
        version, nstr, norb, nelec, total_elec, bond, sum_bond, mode = struct.unpack(
            "<8i", core.read_exact(stream, 32))
        if version != VERSION or mode not in (0, 1):
            raise ValueError("unsupported C++ TNVB request")
        if min(nstr, norb, nelec, bond) <= 0 or nelec % 2 or total_elec < nelec:
            raise ValueError("invalid C++ TNVB dimensions")
        if (total_elec - nelec) % 2 or (mode and sum_bond < bond):
            raise ValueError("invalid inactive electron count or total-state bond")
        structures = read_array(stream, nstr * total_elec, "<i4", (total_elec, nstr))
        overlap = read_array(stream, norb * norb, "<f8", (norb, norb))
        one_electron = read_array(stream, norb * norb, "<f8", (norb, norb))
        npair = norb * (norb + 1) // 2
        eri = read_array(stream, npair * (npair + 1) // 2, "<f8")
        if stream.read(1):
            raise ValueError("trailing C++ TNVB request data")
        if not all(np.all(np.isfinite(v)) for v in (overlap, one_electron, eri)):
            raise ValueError("nonfinite C++ TNVB integrals")
    return (nstr, norb, nelec, total_elec, bond, structures,
            overlap, one_electron, eri), sum_bond, mode


def dense_roots(hamiltonian, overlap):
    core.validate_structure_metric(hamiltonian, overlap)
    norms = np.sqrt(np.diag(overlap))
    h = hamiltonian / norms[:, None] / norms[None, :]
    s = overlap / norms[:, None] / norms[None, :]
    energy, vector = scipy.linalg.eigh(h, s, driver="gvd", check_finite=True)
    coefficients = vector / norms[:, None]
    return energy, coefficients


def evaluate(request, sum_bond, mode):
    # The C++ paired-CSF expansion is invariant under exchanging the two
    # members of an individual singlet pair. The reused XMVB core attaches an
    # oriented-pair phase instead. Canonicalize pair members at this adapter
    # boundary so H/S and reported coefficients follow the C++ convention;
    # preserve pair order, inactive occupations and the caller's raw input.
    nstr, norb, nelec, total_elec, bond, structures, overlap, h1, eri = request
    structures = structures.copy()
    first_active = total_elec - nelec
    for istr in range(nstr):
        active_pairs = structures[first_active:, istr].reshape((-1, 2))
        structures[first_active:, istr] = np.sort(active_pairs, axis=1).ravel()
    request = (nstr, norb, nelec, total_elec, bond, structures, overlap, h1, eri)
    started = time.perf_counter()
    if mode == 0:
        hmat, smat = core.build_hs(request)
        reduced = None
        fit_error = identity_error = rdm_seconds = 0.0
    else:
        result = core.build_hs(request, root=0, sum_bond=sum_bond,
                               d3g_mode=0, d3g_backend=0, sum_mode=0)
        hmat, smat, fit_error, root_energy, root_coefficient, p, pi, d3g, checks = result
        reduced = (p, pi, d3g)
        identity_error = checks["identity_error"]
        rdm_seconds = checks["total_seconds"]
    energy, coefficients = dense_roots(hmat, smat)
    if mode:
        if abs(energy[0] - root_energy) > 1e-10:
            raise RuntimeError("RDM and dense H/S root energies disagree")
        # Preserve exactly the coefficients actually used to construct the RDM.
        coefficients[:, 0] = root_coefficient
    hs_seconds = max(0.0, time.perf_counter() - started - rdm_seconds)
    return (hmat, smat, energy, coefficients, reduced,
            (fit_error, identity_error, hs_seconds, rdm_seconds))


def write_response(path, request, mode, result):
    nstr, norb = request[:2]
    h, s, energy, coefficients, reduced, diagnostics = result
    with open(path, "wb") as stream:
        stream.write(RESPONSE_MAGIC)
        stream.write(struct.pack("<8i", VERSION, 0, nstr, norb, mode, nstr, 0, 0))
        for value in (h, s, energy, coefficients):
            if not np.all(np.isfinite(value)):
                raise RuntimeError("nonfinite C++ TNVB response")
            stream.write(np.asarray(value, dtype="<f8").ravel(order="F").tobytes())
        stream.write(struct.pack("<4d", *diagnostics))
        if mode:
            for value in reduced:
                if not np.all(np.isfinite(value)):
                    raise RuntimeError("nonfinite C++ TNVB density response")
                stream.write(np.asarray(value, dtype="<f8").ravel(order="F").tobytes())


def main(argv):
    if len(argv) != 3:
        raise ValueError("usage: evaluate.py REQUEST RESPONSE")
    request, sum_bond, mode = read_request(argv[1])
    write_response(argv[2], request, mode, evaluate(request, sum_bond, mode))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main(sys.argv))
    except Exception:
        traceback.print_exc()
        raise SystemExit(99)
