# Vendored: `midout` (Gidney, "Inplace Access to the Surface Code Y Basis")

Source: C. Gidney, Data for "Inplace Access to the Surface Code Y Basis", Zenodo,
doi:10.5281/zenodo.7487893 (code.zip, `src/midout/`), licence CC-BY-4.0. Paper: arXiv:2302.07395,
Quantum 8, 1310 (2024).

Copied unchanged except that his `*_test.py` files are left out and one compatibility patch for
numpy >= 2:
`midout/gen/_flow_verifier.py`: `np.bool8` -> `np.bool_` (alias removed in numpy 2).
With this patch and stim 1.16 the generator reproduces the archived noisy circuits byte for byte
(checked by the test item "Gidney's generator reproduces his archived circuits (bridge)" in
`test/test_circuits.jl`).

sha256 of the downloaded archives (2026-11-01):
e72db19c7486a11d0ee588b611d3bfb64c28f05e422eec2dc79e1f272a8f35f9  code.zip
2ded3a58aae29c209a19d539ae7be85659b503f48e899c60449e3336d17596ce  circuits.zip
e35f09cf62cf5fb06966b8e40797780cd8b2f3765684c86fcad6141e0ff4d0bc  stats.csv
