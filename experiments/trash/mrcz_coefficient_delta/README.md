# Archived MRCZ coefficient-delta experiment

This experiment is no longer part of the formal MRCZ implementation or the current paper. Its implementation is archived at `c_codes/trash/coefficient_delta.tpp`; restoring it also requires reconnecting the method declarations, include, CLI routing, and compress/decompress calls removed on 2026-08-29. The result CSV files remain beside the experiment that generated them.

This experiment tests a reversible first-order delta transform across elements after coefficient-major transposition and before Huffman coding. The formal MRCZ bitstream is compared with variants that apply the transform to only the constant mode, all modes with total order at most one, or every mode.

The delta transform operates on quantized `int64_t` coefficients modulo `2^64`, so decompression reconstructs the original coefficient QBs exactly and does not change the error bound. The historical driver is `run_small.py` in this directory.

## Small-data results

The test covers 44 dataset/field/EB combinations: all five fields of `mixlay0`, `T/S00/S01` of `test0`, and `UX/UY/UZ` of `tc_channel0`, each at REL EB `1e-2`, `1e-3`, `1e-4`, and `1e-5`. Every method uses one CPU and one partition. All variants reproduce the formal MRCZ PSNR exactly, and the largest observed `max_error/EB` is `0.997690`.

| Variant | Geometric-mean BR change | Wins | Losses | Best case | Worst case |
|---|---:|---:|---:|---:|---:|
| Constant mode | -7.02% | 31 | 13 | -32.86% | +8.04% |
| Total order at most one | -7.76% | 30 | 14 | -46.52% | +18.27% |
| All modes | +0.63% | 17 | 27 | -49.56% | +61.34% |

The effect is strongly field-dependent. Constant-mode delta has geometric-mean changes of `-7.00%`, `-13.37%`, and `-0.26%` on `mixlay0`, `test0`, and `tc_channel0`; the corresponding order-at-most-one changes are `-10.06%`, `-16.82%`, and `+6.69%`. Unconditional all-mode delta is rejected. A fixed low-order transform is promising but not robust enough for the formal method; any continuation should first predict whether delta is beneficial for each partition or low-order mode using a cheap statistic, then record the selected mask in the bitstream.

`raw_results.csv` contains all BR, PSNR, error, and timing values. `summary.csv` contains the aggregate comparison. A separate 7-partition, 4-CPU `test0/S01` check also reproduced the formal PSNR, MSE, and maximum error exactly.
