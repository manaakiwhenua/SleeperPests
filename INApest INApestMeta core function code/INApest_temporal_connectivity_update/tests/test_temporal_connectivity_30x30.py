#!/usr/bin/env python3
"""Connectivity-level regression tests for the INApest temporal SDD/LDD patch.

These tests use a 30 x 30 (900-node) landscape and reproduce the matrix/array
selection and combination performed by the R functions.  They are deliberately
independent of the model biology so failures isolate temporal-connectivity logic.
"""
from pathlib import Path
import numpy as np

NROW = 30
NCOL = 30
N = NROW * NCOL
T = 4


def validate(sdd, ldd, timesteps):
    if sdd.ndim == 3 and (sdd.shape[0] != sdd.shape[1] or sdd.shape[2] != timesteps):
        raise ValueError("SDDprob 3D array must have dimensions nodes x nodes x Ntimesteps")
    if ldd.ndim == 3 and (ldd.shape[0] != sdd.shape[0] or ldd.shape[1] != sdd.shape[0] or ldd.shape[2] != timesteps):
        raise ValueError("LDDprob 3D array must have dimensions nodes x nodes x Ntimesteps")


def select(x, timestep):
    return x[:, :, timestep] if x.ndim == 3 else x


def combine(sdd, ldd, timestep):
    nsdd = select(sdd, timestep)
    nldd = select(ldd, timestep)
    # The R functions combine matrices as 1-(1-SDD)*(1-LDD).
    return 1.0 - (1.0 - nsdd) * (1.0 - nldd)


def build_landscape():
    sdd = np.zeros((N, N), dtype=np.float32)
    for r in range(NROW):
        for c in range(NCOL):
            i = r * NCOL + c
            if r > 0:
                sdd[i, (r - 1) * NCOL + c] = 0.08
            if r < NROW - 1:
                sdd[i, (r + 1) * NCOL + c] = 0.08
            if c > 0:
                sdd[i, r * NCOL + c - 1] = 0.08
            if c < NCOL - 1:
                sdd[i, r * NCOL + c + 1] = 0.08

    ldd = np.zeros((N, N), dtype=np.float32)
    # One low-probability long-distance edge per node.
    for i in range(N):
        j = (i + N // 2 + 17) % N
        if j != i:
            ldd[i, j] = 0.015
    return sdd, ldd


def assert_raises(fn):
    try:
        fn()
    except ValueError:
        return
    raise AssertionError("Expected ValueError")


def main():
    sdd, ldd = build_landscape()
    sdd_repeat = np.repeat(sdd[:, :, None], T, axis=2)
    ldd_repeat = np.repeat(ldd[:, :, None], T, axis=2)

    # Dynamic copies: timestep 2 (R timestep index 2; Python index 1) closes the
    # central vertical SDD connection. Timestep 3 changes LDD connectivity.
    sdd_dynamic = sdd_repeat.copy()
    for r in range(NROW):
        left = r * NCOL + (NCOL // 2 - 1)
        right = r * NCOL + (NCOL // 2)
        sdd_dynamic[left, right, 1] = 0.0
        sdd_dynamic[right, left, 1] = 0.0

    ldd_dynamic = ldd_repeat.copy()
    ldd_dynamic[:, :, 2] = 0.0
    for i in range(N):
        j = (i + N // 3 + 11) % N
        if j != i:
            ldd_dynamic[i, j, 2] = 0.025

    # Shape validation.
    validate(sdd, ldd, T)
    validate(sdd_repeat, ldd_repeat, T)
    validate(sdd_dynamic, ldd, T)
    validate(sdd, ldd_dynamic, T)
    assert_raises(lambda: validate(np.zeros((N, N - 1, T), dtype=np.float32), ldd_repeat, T))
    assert_raises(lambda: validate(np.zeros((N, N, T - 1), dtype=np.float32), ldd_repeat, T))
    assert_raises(lambda: validate(sdd, np.zeros((N, N, T - 1), dtype=np.float32), T))

    # Backward compatibility: repeating a static matrix in every 3D slice must
    # produce exactly the same selected matrices and combined probabilities.
    static_combined = combine(sdd, ldd, 0)
    for t in range(T):
        np.testing.assert_array_equal(select(sdd_repeat, t), sdd)
        np.testing.assert_array_equal(select(ldd_repeat, t), ldd)
        np.testing.assert_array_equal(combine(sdd_repeat, ldd_repeat, t), static_combined)

    # Mixed 2D/3D inputs must select only the temporal argument.
    for t in range(T):
        np.testing.assert_array_equal(combine(sdd_repeat, ldd, t), static_combined)
        np.testing.assert_array_equal(combine(sdd, ldd_repeat, t), static_combined)

    # Temporal SDD should affect only the intended timestep.
    np.testing.assert_array_equal(combine(sdd_dynamic, ldd, 0), static_combined)
    assert not np.array_equal(combine(sdd_dynamic, ldd, 1), static_combined)
    np.testing.assert_array_equal(combine(sdd_dynamic, ldd, 2), static_combined)
    np.testing.assert_array_equal(combine(sdd_dynamic, ldd, 3), static_combined)

    # Temporal LDD should likewise affect only the intended timestep.
    np.testing.assert_array_equal(combine(sdd, ldd_dynamic, 0), static_combined)
    np.testing.assert_array_equal(combine(sdd, ldd_dynamic, 1), static_combined)
    assert not np.array_equal(combine(sdd, ldd_dynamic, 2), static_combined)
    np.testing.assert_array_equal(combine(sdd, ldd_dynamic, 3), static_combined)

    # Direct functional check across the closed central edge.
    source = 10 * NCOL + (NCOL // 2 - 1)
    target = 10 * NCOL + (NCOL // 2)
    assert static_combined[source, target] > 0
    assert combine(sdd_dynamic, ldd, 1)[source, target] == 0

    # Establishment weighting follows the same target-column sweep as R.
    env = np.linspace(0.4, 1.0, N, dtype=np.float32)
    bpam_static = static_combined * env[None, :]
    bpam_repeat = combine(sdd_repeat, ldd_repeat, 3) * env[None, :]
    np.testing.assert_array_equal(bpam_static, bpam_repeat)
    bpam_barrier = combine(sdd_dynamic, ldd, 1) * env[None, :]
    assert bpam_barrier[source, target] == 0

    # All probabilities remain finite and within [0,1].
    for m in (static_combined, combine(sdd_dynamic, ldd_dynamic, 1), combine(sdd_dynamic, ldd_dynamic, 2)):
        assert np.isfinite(m).all()
        assert float(m.min()) >= 0.0
        assert float(m.max()) <= 1.0

    # Report useful landscape-level diagnostics.
    static_edges = int(np.count_nonzero(static_combined))
    barrier_edges = int(np.count_nonzero(combine(sdd_dynamic, ldd, 1)))
    ldd_changed_edges = int(np.count_nonzero(combine(sdd, ldd_dynamic, 2)))
    print("PASS: 30x30 temporal-connectivity regression tests")
    print(f"nodes={N}, timesteps={T}")
    print(f"static combined non-zero links={static_edges}")
    print(f"barrier timestep non-zero links={barrier_edges}")
    print(f"changed-LDD timestep non-zero links={ldd_changed_edges}")
    print("verified: 2D == repeated 3D, mixed 2D/3D, SDD temporal change, LDD temporal change, validation, establishment weighting")


if __name__ == "__main__":
    main()
