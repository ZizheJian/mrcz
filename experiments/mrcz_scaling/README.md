# MRCZ scaling experiments

## Formal weak scaling after the large-data RD experiment

`run_weak_scaling.py` uses the completed full-data RD curves to compare formal MRCZ and FDWZ at fixed quality. pb0 targets 100 dB and uses 1--64 CPUs; jicf targets 150 dB and uses 1--128 CPUs. Every method/field has its own REL EB interpolated in PSNR--log10(REL EB) space and a corresponding ABS EB computed from the full-field range. The range tests use that fixed ABS EB so selecting a local element interval cannot silently change the effective error bound.

The maximum CPU count processes the complete dataset with one global partition per CPU. Lower CPU counts process groups from the same fixed global partition grid. Each CPU point aggregates four observations: one available group is repeated four times, two groups are each repeated twice, and four or more available groups use four spatially distributed groups. Across observations and fields, throughput is total original bytes divided by total internal backend time; distinct physical regions are not summarized with a median. Weighted MSE is converted to PSNR with the full-field range.

- `run_weak_pb0.sbatch`: pb0, 64 CPUs maximum, 100 dB.
- `run_weak_jicf.sbatch`: jicf, 128 CPUs maximum, 150 dB.
- `weak_<dataset>_raw.csv`: every method/field/CPU/global-partition-group observation.
- `weak_<dataset>_summary.csv`: aggregate compression/decompression throughput and weak-scaling efficiency.
- `weak_<dataset>_quality.csv`: aggregate BR, MSE, PSNR, maximum error, and calibrated EB for every field.

## Preliminary combined strong/weak experiment

This experiment compares the current formal `MRCZ-3D-GLL` and `SZ3-3D-GLL-optimized` implementations on `jicf0/S00`.

- Strong scaling uses REL EB `1e-6`, elements `[1350000,1862000)`, and a fixed 128-partition bitstream for every CPU count.
- Weak scaling uses fixed ABS EB `5e-5`, 20000 elements and one partition per CPU, beginning at element 300000. ABS mode prevents the changing dynamic range of the growing element range from changing the effective error bound.
- CPU counts are 1, 2, 4, 8, 16, 32, 64, and 128. Each case is repeated three times in interleaved order.
- `OMP_PLACES=cores` and `OMP_PROC_BIND=close` make placement reproducible. Phase profiling is disabled.
- `raw_results.csv` stores every repetition. `summary.csv` stores medians, aggregate throughput, speedup, and efficiency. All throughput values come from the compressor's internal backend timers.

For strong scaling, speedup and efficiency are relative to the one-CPU measurement. For weak scaling, the reported efficiency is the retained per-CPU throughput relative to one CPU.

The first submitted run is Slurm job `20143133`.
