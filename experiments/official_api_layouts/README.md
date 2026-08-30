# Official API tensorization experiment

This experiment tests direct official compressor APIs plus MRCZ and the optimized GLL predictor on three small NekRS datasets, two MFEM examples, and one Nektar++ example. Every method uses one CPU and one partition covering the complete selected field.

Run `generate_mfem_and_nektar_data.py` from the project root when the MFEM and Nektar++ validation inputs need to be regenerated.

The tested layouts are:

- SZ3, QoZ, and ZFP: flattened 1D, `[element,N]` 2D, `[element*N0,N1,N2]` fused 3D, one 3D API call per element, and `[element,N0,N1,N2]` 4D.
- SZ2: the same five API calls. The official default regression implementation internally merges the first two dimensions of 4D data, so its fused-3D and 4D results may be identical.
- SPERR: `[element,N^3]` 2D, fused 3D, and one 3D API call per element. Its official C API does not provide 1D or 4D compression.
- QoZ removes singleton dimensions before calling the API. Its q4 tuning is only active for 2D and 3D; the 1D layout uses the official interpolation backend because the current official `ALGO_INTERP_LORENZO` path corrupts the heap in 1D.

Per-element methods store one 64-bit stream length per element plus an element count, and this necessary metadata is included in BR. `run_small_rd.py` records REL EB values from `1e-2` through `1e-6` in `<data_path>_<field>.csv`, records backend throughput in `<data_path>_<field>_throughput.csv`, and supports resuming completed dataset/method/EB combinations from those files.

Run `python3 run_small_rd.py` after building and installing MRCZ, then run `python3 summarize_rd.py`. The summary reads the per-field CSV files, interpolates BR inside the PSNR range shared by every layout of one compressor, and reports each layout's BR regret relative to the pointwise-best layout. The known non-discriminating `test0/UZ` case remains in its per-field CSV but is excluded from layout selection.

The small-data layouts selected for the formal large-data experiment are:

- SZ3: 4D.
- SZ2: fused 3D. The official library removes singleton dimensions first. For full 3D elements, its 4D path merges `E` and `N0` internally and is bitstream-identical to fused 3D, so only fused 3D is retained.
- QoZ: 4D.
- ZFP: fused 3D.
- SPERR: 2D.

SZ3's official float path showed a maximum numerical overshoot of about 0.031% when given the requested EB directly, so its wrapper passes `0.999*abs_eb`; the regenerated results all satisfy the requested bound. Official SZ2's per-element 3D path aborts inside its allocator for the MFEM ex1 `5x5x1` elements at REL `1e-4`; this failure remains documented here and in `raw.log`.
