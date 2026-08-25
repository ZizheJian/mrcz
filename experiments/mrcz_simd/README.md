# MRCZ SIMD experiments

`run_profile_large.sbatch` profiles the current MRCZ implementation on one 100000-element `pb0` partition. `MRCZ_PROFILE_PHASES=1` enables program-internal accumulated partition CPU-time reports for transforms, quantization, Huffman, and Zstd; normal runs remain unchanged when the variable is absent. Jobs `20141358` and `20141422` captured the scalar and pointer-refactored versions respectively.

The GCC 11.2 vectorization report generated with the release flags (`-Ofast -march=native`) reported zero vectorized loops in the MRCZ transform and quantization path. Tensor transforms were rejected mainly because of nested accumulation loops and gather/strided loads. Quantization was rejected mainly around floating-point rounding and conversion to `int64_t`.

## 2026-08-25 baseline profile

`pb0/T`, elements `[1000000,1100000)`, REL `1e-4`, one 100000-element partition, one application CPU:

- Compression: 2.896769 seconds. Forward transform 27.66%, inverse transform 27.86%, coefficient quantization 3.24%, residual quantization 2.78%, coefficient Huffman 31.82%, residual Huffman 5.03%, Zstd 1.62%.
- Decompression: 1.991796 seconds. Inverse transform 44.43%, coefficient Huffman 45.93%, residual Huffman 4.40%, coefficient reconstruction 1.87%, residual reconstruction 2.87%, Zstd 0.51%.

A first `#pragma omp simd` experiment annotated the six contiguous innermost transform loops without changing their order. On `pb0/T`, 10000 elements, one CPU, the three scalar baseline compression times were 0.310675, 0.317815, and 0.316758 seconds; the annotated times were 0.324077, 0.321524, and 0.321048 seconds. The compiler still did not report the transform loops as successfully vectorized. The annotations were removed because this approach provided no speedup.

Replacing repeated `vector[]` index expressions with precomputed contiguous input/output pointers and scalar weights enabled GCC to vectorize all six transform inner loops with 32-byte vectors. This preserves the arithmetic loop order and does not use compiler-specific intrinsics. On the same 10000-element login-node sample, three compression times were 0.176330, 0.176358, and 0.176711 seconds, about 1.78x faster than the scalar baseline.

The 100000-element compute-node rerun produced exactly the same BR, PSNR, and maximum error as the baseline. Compression time fell from 2.896769 to 2.126830 seconds (26.57% less time; throughput 67.424261 to 91.832669 MB/s, 1.362x). Decompression time fell from 1.991796 to 1.440229 seconds (27.69% less time; throughput 98.058481 to 135.612136 MB/s, 1.383x). Forward transform fell from 0.773426 to 0.372865 seconds (2.074x), compression-side inverse transform from 0.779152 to 0.391879 seconds (1.988x), and decompression-side inverse transform from 0.879212 to 0.365360 seconds (2.406x).

The four-CPU `jicf/S00` job `20141452` completed all nine range/error-bound combinations in 38 seconds with exit code 0. It confirms that the vectorized implementation works with partition-level OpenMP. Its endpoint times must not be compared directly with the older coefficient-tree log because this job enabled per-element phase timers while the older scalar job did not. The phase results nevertheless show that coefficient Huffman accounts for 52%--78% of measured work in the tight-EB cases, so additional handwritten transform AVX2 would have a limited endpoint ceiling there. If a formal multicore speedup number is needed, scalar and pointer-refactored implementations should be interleaved in the same allocation with phase profiling disabled.

Decision: keep the pointer refactor and compiler-generated AVX2, and do not add handwritten AVX2 intrinsics to the transform. Reconsider intrinsics only if a later assembly/profile audit finds an important transform loop that is no longer vectorized on a target compiler.

## Cross-element batched prototype

`MRCZ-3D-GLL-element-batched` preserves the formal `MRCZ-3D-GLL` path and adds a separate experiment. It packs the same local position from eight float elements or four double elements into one SIMD lane group, applies all three forward/inverse transform stages across element lanes, and directly stores coefficient QBs in the coefficient-major order consumed by Huffman. The latter removes the formal path's full-partition coefficient transpose. Incomplete batches use the formal scalar-element transform.

GCC reports all six batched transform lane loops as vectorized with 32-byte vectors. On `test0/T`, including a separate 509-element/three-partition tail test, the field bitstream, container, BR, PSNR, MSE, and maximum error are byte-for-byte or numerically identical to `MRCZ-3D-GLL`. The full 512-element one-CPU profile reduced compression time from 0.004857 to 0.003190 seconds and decompression from 0.002315 to 0.001555 seconds, but these endpoint ratios favor batching because phase profiling calls the clock once per element in the formal path and once per batch in the prototype.

On double-precision `tc_channel0/UX`, bitstreams and errors were also identical. Batched forward/inverse transforms were slightly slower than the formal element-internal SIMD transforms, while removing coefficient transpose reduced the lossless coefficient stage enough to change compression from 0.009046 to 0.008677 seconds and decompression from 0.006238 to 0.004442 seconds.

With phase profiling disabled, five interleaved login-node `pb0/T` runs on 10000 elements had median compression times of 0.176441 seconds for `MRCZ-3D-GLL` and 0.153003 seconds for the batched prototype, a 13.29% time reduction or 1.153x speedup. The node was noisy, but four of five paired runs favored batching by 8%--15%; the fifth was an outlier with a larger difference. Stage profiles indicate that the transform itself is not faster on this 8-cubed shape; the observed endpoint gain primarily comes from avoiding the naive full-partition coefficient transpose.

Large job `20141698` did not provide a valid comparison. On node `a120`, the batched compressor entered Linux `D` state while reading the shared filesystem; `sstat` showed only about 19 CPU seconds after nearly five wall-clock minutes, and process inspection confirmed uninterruptible I/O wait. The preceding formal-path run was also abnormally slow (19.09 seconds versus the usual 2--3 seconds). Its absolute and relative timings must not be used. A future large comparison should stage a prepared subset on node-local storage or otherwise avoid repeated reads from the congested shared file.

## Cross-element/transpose ablation

Two additional experimental methods isolate the two changes that were originally combined in `MRCZ-3D-GLL-element-batched`:

- `MRCZ-3D-GLL-element-batched-original-transpose` uses the cross-element transform but writes element-major coefficient QBs and retains the original full-partition transpose.
- `MRCZ-3D-GLL-tiled-transpose` retains the formal element-internal transform and replaces only the naive transpose with a 32-by-32 cache-blocked transpose.

All four methods produced byte-identical field bitstreams and containers on float `test0/T` with 509 elements and three partitions, including incomplete SIMD batches. They also produced byte-identical outputs on double `tc_channel0/UX` under the same tail and partition test.

On `pb0/T`, elements `[1000000,1010000)`, REL `1e-4`, one partition and one CPU, seven interleaved runs with phase profiling disabled gave these median backend times:

| Method | Compression (s) | Decompression (s) | Compression change | Decompression change |
| --- | ---: | ---: | ---: | ---: |
| Formal element-internal SIMD + naive transpose | 0.170182 | 0.116888 | baseline | baseline |
| Cross-element SIMD + naive transpose | 0.184100 | 0.120733 | +8.18% | +3.29% |
| Element-internal SIMD + tiled transpose | 0.150987 | 0.099728 | -11.28% | -14.68% |
| Cross-element SIMD + direct coefficient-major output | 0.148497 | 0.101906 | -12.74% | -12.82% |

Three profiled runs explain the result. The original coefficient lossless stage took median 0.062851 seconds during compression and 0.062105 seconds during decompression. Tiling reduced these to 0.041073 and 0.045914 seconds. Direct coefficient-major output reduced them further to 0.025421 and 0.025572 seconds, but cross-element forward/inverse transforms were slower and additionally required element packing. For example, compression-side forward/inverse medians changed from 0.032280/0.032185 seconds to 0.039403/0.039270 seconds before packing.

Conclusion: vectorizing across elements is not beneficial for the current 8-cubed MRCZ transform. Cache-blocking the coefficient transpose provides nearly the same combined compression-plus-decompression benefit with much less code and no packing. The 32-by-32 tiled transpose has therefore been promoted to the formal `MRCZ-3D-GLL` path, and the redundant `MRCZ-3D-GLL-tiled-transpose` method name has been removed. The cross-element implementation remains only as an experimental record. Larger and multicore validation belongs to the formal throughput/scaling experiment rather than further SIMD development.
