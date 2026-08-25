# MRCZ coefficient multi-tree experiments

These experiments tested whether coefficient modes with different total polynomial orders should use separate Huffman trees. All trees in one partition were concatenated before a shared Zstd compression step, so reported BR includes every Huffman tree's metadata.

## Methods

- Single tree: formal `MRCZ-3D-GLL`, with coefficient-major ordering.
- Two trees: modes with total order at most a threshold versus the remaining modes; thresholds 0 through 4 were tested.
- Three trees: constant, total-order-1 linear, and remaining modes.
- Four trees: constant, linear, total-order-2 quadratic, and remaining modes.

## Small-data result

On test0/T, mixlay/T, and tc_channel/UX with one CPU and REL EB from `1e-2` through `1e-5`, three trees reduced BR by a geometric mean of 8.54% relative to one tree. Four trees reduced BR by 9.67%, but the additional improvement over three trees was only 1.24% and was not stable across all cases. PSNR was identical because only lossless coefficient coding changed.

## Large-partition result

`run_large_partition.sbatch` tested jicf/S00 on three separated ranges of 100000 elements, using four CPUs and therefore approximately 25000 elements per partition. REL EB values were `1e-4`, `1e-6`, and `1e-8`.

- Three trees versus one tree: geometric-mean BR change `-1.13%`.
- Four trees versus one tree: geometric-mean BR change `-1.14%`.
- Four trees versus three trees: geometric-mean BR change `-0.012%`.
- The three sampled ranges gave three-tree BR changes of `-0.76%`, `+0.15%`, and `-2.77%` relative to one tree.

The large-partition benefit was small and region-dependent, while timing did not show a stable trend. Formal MRCZ therefore retains one coefficient Huffman tree. The implementation is archived in `c_codes/trash/coefficient_multitree_experiments.tpp`.

The exact large-partition output is in `slurm-20140938.out`; `slurm-20140938.err` is empty.
