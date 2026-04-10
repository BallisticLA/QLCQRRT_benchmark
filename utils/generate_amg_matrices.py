#!/usr/bin/env python3
"""Generate coupled (K, V) matrix pairs from PyAMG for generalized SVD/LS benchmarks.

K = fine-grid stiffness matrix (m x m, sparse SPD, not diagonal)
V = prolongation operator  (m x n, tall sparse, coupled to K)

Both matrices arise from the same algebraic multigrid hierarchy, so they are
inherently coupled by the discretization.

Column scaling of V controls the condition number of L^{-1}V (where LL^T = K):
  V_scaled = V @ diag(sigma),  sigma = logspace(0, log10(kappa_mult), n)
This corresponds to rescaling the coarse-space basis functions, which is
physically meaningful and preserves the coupled structure.

Column permutation (--permute) randomly reorders V's columns before scaling.
This distributes scale factors across the aggregate structure rather than
aligning them monotonically, breaking the near-orthogonality that AMG
prolongation columns naturally have. This is needed to create examples
where CholQR fails.

Usage:
    python generate_amg_matrices.py [--grid SIZES] [--kappa MULTS] [--permute] [--seed SEED]

    --grid    Space-separated grid sizes (default: 32 64 128)
    --kappa   Space-separated condition number multipliers (default: 1 1e3 1e6 1e9)
              1 = unscaled (well-conditioned baseline)
    --permute Randomly permute columns before scaling
    --seed    Random seed for permutation (default: 42)

Examples:
    python generate_amg_matrices.py
    python generate_amg_matrices.py --grid 64 --kappa 1 1e6 1e9
    python generate_amg_matrices.py --grid 64 --kappa 1e4 1e6 1e8 --permute

Output goes to ./input_matrices/ (relative to this script's directory).
"""

import pyamg
import pyamg.gallery
import scipy.io
import scipy.sparse
import scipy.linalg
import numpy as np
import os
import sys
import argparse


def generate_poisson_problem(grid_size, kappa_mult, output_dir, permute=False, seed=42):
    """Generate 2D Poisson problem at given grid size with optional column scaling.

    Args:
        grid_size: Number of grid points per dimension (total DOFs = grid_size^2).
        kappa_mult: Condition number multiplier for column scaling.
                    1.0 = unscaled (baseline).
                    >1   = columns scaled by logspace(0, log10(kappa_mult), n).
        output_dir: Directory to write .mtx files.
        permute: If True, randomly permute columns before scaling.
        seed: Random seed for permutation.

    Returns:
        dict with keys: grid, kappa_mult, m, n, nnz_K, nnz_V, K_file, V_file, kappa_LiV.
    """
    # 2D Poisson with finite-element stencil
    stencil = pyamg.gallery.diffusion_stencil_2d(type='FE', epsilon=1.0)
    K = pyamg.gallery.stencil_grid(stencil, (grid_size, grid_size), format='csr')
    m = K.shape[0]

    # Build smoothed aggregation AMG hierarchy
    ml = pyamg.smoothed_aggregation_solver(K, max_coarse=10)

    # V = prolongation operator from coarse level 1 to fine level 0
    V = ml.levels[0].P
    n = V.shape[1]

    # Optionally permute columns before scaling
    if permute and kappa_mult > 1.0:
        rng = np.random.default_rng(seed)
        perm = rng.permutation(n)
        V = V.tocsc()[:, perm].tocsr()

    # Apply column scaling if kappa_mult > 1
    if kappa_mult > 1.0:
        sigma = np.logspace(0, np.log10(kappa_mult), n)
        D = scipy.sparse.diags(sigma)
        V = V @ D

    # Verify properties
    assert K.shape == (m, m), f"K must be square, got {K.shape}"
    assert V.shape == (m, n), f"V shape mismatch: {V.shape}"
    assert n < m, f"V must be tall: m={m}, n={n}"

    # Check K is SPD (quick check: all diagonal entries positive)
    diag = K.diagonal()
    assert np.all(diag > 0), "K diagonal has non-positive entries"

    # Compute actual κ(L^{-1}V) for small enough matrices
    kappa_LiV = None
    if m <= 20000:
        K_dense = K.toarray()
        V_dense = V.toarray()
        L = scipy.linalg.cholesky(K_dense, lower=True)
        LiV = scipy.linalg.solve_triangular(L, V_dense, lower=True)
        svals = np.linalg.svd(LiV, compute_uv=False)
        kappa_LiV = svals[0] / svals[-1]

    # File naming: include kappa_mult and permute in name
    if kappa_mult <= 1.0:
        kappa_tag = "base"
    else:
        perm_tag = "p" if permute else ""
        kappa_tag = f"{perm_tag}k{kappa_mult:.0e}"
    prefix = f"poisson2d_{grid_size}x{grid_size}_{kappa_tag}"
    K_file = os.path.join(output_dir, f"{prefix}_K.mtx")
    V_file = os.path.join(output_dir, f"{prefix}_V.mtx")

    # Write K as symmetric (only lower triangle stored)
    scipy.io.mmwrite(K_file, scipy.sparse.tril(K), symmetry='symmetric')
    scipy.io.mmwrite(V_file, V)

    perm_str = " (permuted)" if permute and kappa_mult > 1.0 else ""
    kappa_str = f"{kappa_LiV:.2e}" if kappa_LiV else "not computed"
    print(f"  Poisson {grid_size}x{grid_size}, kappa_mult={kappa_mult:.0e}{perm_str}:")
    print(f"    K: {m} x {m}, nnz = {K.nnz}")
    print(f"    V: {m} x {n}, nnz = {V.nnz}")
    print(f"    κ(L⁻¹V) = {kappa_str}")
    print(f"    Files: {os.path.basename(K_file)}, {os.path.basename(V_file)}")

    return {
        'grid': grid_size, 'kappa_mult': kappa_mult, 'permute': permute,
        'm': m, 'n': n,
        'nnz_K': K.nnz, 'nnz_V': V.nnz,
        'K_file': os.path.basename(K_file),
        'V_file': os.path.basename(V_file),
        'kappa_LiV': kappa_LiV,
    }


def generate_random_problem(grid_size, kappa_mult, output_dir, seed=42):
    """Generate 2D Poisson K with a dense random V having controlled conditioning.

    V = U_rand @ diag(sigma) @ Vt_rand where U_rand, Vt_rand are random orthogonal.
    This produces a dense V with κ(V) ≈ kappa_mult and κ(L^{-1}V) ≈ kappa_mult.

    At κ(L^{-1}V) ≈ 1e7, CholQR fails (~1e-3 orthogonality error) while
    CQRRT_linop (~1e-9) and sCholQR3 (~1e-10) succeed.

    Args:
        grid_size: Grid points per dimension (total DOFs = grid_size^2).
        kappa_mult: Target condition number for V's singular values.
        output_dir: Directory to write .mtx files.
        seed: Random seed for orthogonal factors.

    Returns:
        dict with matrix info.
    """
    stencil = pyamg.gallery.diffusion_stencil_2d(type='FE', epsilon=1.0)
    K = pyamg.gallery.stencil_grid(stencil, (grid_size, grid_size), format='csr')
    m = K.shape[0]

    # Use AMG hierarchy to determine n (number of columns)
    ml = pyamg.smoothed_aggregation_solver(K, max_coarse=10)
    n = ml.levels[0].P.shape[1]

    rng = np.random.default_rng(seed)
    U_rand, _ = np.linalg.qr(rng.standard_normal((m, n)))
    Vt_rand, _ = np.linalg.qr(rng.standard_normal((n, n)))

    if kappa_mult <= 1.0:
        sigma_v = np.ones(n)
        kappa_tag = "rbase"
    else:
        sigma_v = np.logspace(0, np.log10(kappa_mult), n)
        kappa_tag = f"rk{kappa_mult:.0e}"

    V_dense = U_rand @ np.diag(sigma_v) @ Vt_rand
    V_sparse = scipy.sparse.csr_matrix(V_dense)

    # Compute actual κ(L^{-1}V) for feasible sizes
    kappa_LiV = None
    if m <= 20000:
        K_dense = K.toarray()
        L = scipy.linalg.cholesky(K_dense, lower=True)
        LiV = scipy.linalg.solve_triangular(L, V_dense, lower=True)
        svals = np.linalg.svd(LiV, compute_uv=False)
        kappa_LiV = svals[0] / svals[-1]

    prefix = f"random_{grid_size}x{grid_size}_{kappa_tag}"
    K_file = os.path.join(output_dir, f"{prefix}_K.mtx")
    V_file = os.path.join(output_dir, f"{prefix}_V.mtx")

    scipy.io.mmwrite(K_file, scipy.sparse.tril(K), symmetry='symmetric')
    scipy.io.mmwrite(V_file, V_sparse)

    nnz_V = V_sparse.nnz
    kappa_str = f"{kappa_LiV:.2e}" if kappa_LiV else "not computed"
    print(f"  Random {grid_size}x{grid_size}, kappa_mult={kappa_mult:.0e}:")
    print(f"    K: {m} x {m}, nnz = {K.nnz}")
    print(f"    V: {m} x {n}, nnz = {nnz_V} (dense)")
    print(f"    κ(L⁻¹V) = {kappa_str}")
    print(f"    Files: {os.path.basename(K_file)}, {os.path.basename(V_file)}")

    return {
        'grid': grid_size, 'kappa_mult': kappa_mult, 'permute': False,
        'm': m, 'n': n,
        'nnz_K': K.nnz, 'nnz_V': nnz_V,
        'K_file': os.path.basename(K_file),
        'V_file': os.path.basename(V_file),
        'kappa_LiV': kappa_LiV,
    }


def write_metadata(output_dir, entries):
    """Write metadata file listing all generated matrix pairs."""
    meta_file = os.path.join(output_dir, "metadata.txt")
    with open(meta_file, 'w') as f:
        f.write(f"# PyAMG-generated matrix pairs for GSVD/LS benchmark\n")
        f.write(f"# num_pairs\n")
        f.write(f"{len(entries)}\n")
        f.write(f"# grid  kappa_mult  m  n  nnz_K  nnz_V  kappa_LiV  K_file  V_file\n")
        for e in entries:
            kv = f"{e['kappa_LiV']:.6e}" if e['kappa_LiV'] else "NA"
            f.write(f"{e['grid']}  {e['kappa_mult']:.0e}  {e['m']}  {e['n']}  "
                    f"{e['nnz_K']}  {e['nnz_V']}  {kv}  "
                    f"{e['K_file']}  {e['V_file']}\n")
    print(f"\nMetadata written to {meta_file}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Generate coupled (K, V) matrix pairs from PyAMG.')
    parser.add_argument('--grid', nargs='+', type=int, default=[32, 64, 128],
                        help='Grid sizes (default: 32 64 128)')
    parser.add_argument('--kappa', nargs='+', type=float, default=[1, 1e3, 1e6, 1e9],
                        help='Condition number multipliers (default: 1 1e3 1e6 1e9)')
    parser.add_argument('--mode', choices=['amg', 'random'], default='amg',
                        help='Matrix generation mode: amg (sparse prolongation) or random (dense)')
    parser.add_argument('--permute', action='store_true',
                        help='[amg mode] Randomly permute columns before scaling')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for permutation/random generation (default: 42)')

    # Also support legacy positional args: generate_amg_matrices.py 32 64 128
    if len(sys.argv) > 1 and sys.argv[1].isdigit():
        grid_sizes = [int(x) for x in sys.argv[1:]]
        kappa_mults = [1.0]
        mode = 'amg'
        permute = False
        seed = 42
    else:
        args = parser.parse_args()
        grid_sizes = args.grid
        kappa_mults = args.kappa
        mode = args.mode
        permute = args.permute
        seed = args.seed

    # Output directory: input_matrices/ next to this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, "input_matrices")
    os.makedirs(output_dir, exist_ok=True)

    print(f"Mode: {mode}")
    print(f"Grid sizes: {grid_sizes}")
    print(f"Kappa multipliers: {kappa_mults}")
    print(f"Output directory: {output_dir}\n")

    entries = []
    for gs in grid_sizes:
        for km in kappa_mults:
            if mode == 'random':
                entry = generate_random_problem(gs, km, output_dir, seed=seed)
            else:
                entry = generate_poisson_problem(gs, km, output_dir,
                                                 permute=permute, seed=seed)
            entries.append(entry)
            print()

    write_metadata(output_dir, entries)
    print("Done.")
