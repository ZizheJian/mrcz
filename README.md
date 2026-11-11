# MRCZ: Error-Bounded Lossy Compression for High-Order Element Data

**English** | [中文](README_CN.md)

This repository provides two error-bounded lossy compressors for data produced by high-order finite and spectral element simulations, such as NekRS/Nek5000, MFEM, and Nektar++. It also contains the benchmark driver, baseline integrations, experiment configurations, and plotting scripts used in our paper.

High-order element data consist of many elements with irregular global connectivity. Within each element, values are sampled on a small tensor-product grid of Gauss--Lobatto--Legendre (GLL) nodes that is shared by all elements. The two compressors exploit this organization:

- **FDWZ** (Fast Distance-Weighted Interpolation-Based Lossy Compressor) predicts each element with multilevel interpolation weighted by GLL node distances and encodes different prediction levels with separate Huffman trees. It targets high throughput.
- **MRCZ** (Model-Representation-Based Compressor) represents each element with transform coefficients, corrects the transform prediction with pointwise residuals, and groups corresponding coefficients across elements for more effective lossless coding. It targets higher compression ratios.

## Features

- Strict pointwise error bounds in absolute (`ABS`) or value-range-relative (`REL`) mode.
- Native readers for NekRS, MFEM, and Nektar++ output.
- Partition-level OpenMP parallelism, where each partition is compressed into an independent bitstream.
- Command-line tools and a header-only C++ core that can be imported through CMake.
- Unmodified official APIs of SZ3, SZ2, QoZ 2.0 (HPEZ), ZFP, SPERR, and MGARD-X as baselines, with six candidate data layouts.

## Repository Layout

| Path | Contents |
|---|---|
| `c_codes/` | C++ source of FDWZ, MRCZ, baseline wrappers, and the command-line tools |
| `configs/` | JSON experiment configurations |
| `py_funcs/` | Analysis, summary, and plotting scripts |
| `results/` | Experiment CSV files and intermediate plotting data |
| `third_party/` | Git submodules for baselines, simulation frameworks, and data readers |
| `t1_start_c_code.py` | Main experiment driver |
| `generate_mfem_and_nektar_data.py` | Generates the MFEM and Nektar++ test data |
| `run_*.sbatch` | Slurm scripts used for the large-scale experiments |

## Requirements

- Python 3, CMake, and a C++ compiler with C++20 and OpenMP support.
- NumPy, SciPy, Matplotlib, Pillow, and Plotly for analysis and plotting.
- Protobuf for building MGARD-X.

## Building

Clone the repository with its submodules:

```bash
git clone --recurse-submodules git@github.com:ZizheJian/mrcz.git
cd mrcz
```

For an existing clone, initialize the submodules:

```bash
git submodule update --init --recursive
```

Create the Python environment:

```bash
python3 -m venv mrcz_env
source mrcz_env/bin/activate
python -m pip install --upgrade pip
python -m pip install numpy scipy matplotlib pillow plotly
```

MGARD-X requires Protobuf. On Purdue Anvil, load it before building:

```bash
module load protobuf/3.11.4
```

The experiment driver builds all dependencies automatically. To build them manually:

```bash
cmake -S third_party/sperr -B build/sperr_build -DCMAKE_INSTALL_PREFIX="$PWD/install/sperr_install" -DBUILD_SHARED_LIBS=OFF
cmake --build build/sperr_build --parallel 8
cmake --install build/sperr_build

cmake -S third_party/mfem -B build/mfem_build -DCMAKE_INSTALL_PREFIX="$PWD/install/mfem_install"
cmake --build build/mfem_build --parallel 8
cmake --install build/mfem_build

cmake -S third_party/MGARD -B build/mgard_serial -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/install/mgard_install" -DMGARD_ENABLE_SERIAL=ON -DMGARD_ENABLE_OPENMP=OFF -DMGARD_ENABLE_CUDA=OFF -DMGARD_ENABLE_HIP=OFF -DMGARD_ENABLE_SYCL=OFF -DMGARD_ENABLE_MDR=OFF -DMGARD_ENABLE_UNSTRUCTURED=OFF -DMGARD_ENABLE_BENCHMARKS=OFF -DMGARD_ENABLE_DOCS=OFF -DBUILD_TESTING=OFF
cmake --build build/mgard_serial --parallel 8
cmake --install build/mgard_serial

cmake -S c_codes -B build -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build --parallel 8
cmake --install build
```

The command-line tools are installed as `install/bin/mrcz_compress` and `install/bin/mrcz_decompress`.

## Quick Start

### Command-Line Tools

The following example compresses and decompresses two fields of a small NekRS dataset included in a submodule, using FDWZ with a relative error bound of 1e-4 on 8 threads:

```bash
install/bin/mrcz_compress -format nekrs -i third_party/pySEMTools/examples/data/mixlay0.f00001 -fields '[UX,P]' -M REL 1e-4 -method FDWZ -max-cpus 8
```

```bash
install/bin/mrcz_decompress -format nekrs -i third_party/pySEMTools/examples/data/mixlay0.f00001 -fields '[UX,P]' -M REL 1e-4 -method FDWZ -max-cpus 8
```

Each compressed field is written to `<input>_<field>.<method>`, together with the metadata files `<input>_header.mrcz` and `<input>_element_map.mrcz`. The decompressor writes the reconstructed data to `<input>.<method>.bin` and reports the compression ratio, bit rate, maximum error, and PSNR. Both tools accept the same options:

| Option | Meaning |
|---|---|
| `-format nekrs\|mfem\|nektar++` | Input format |
| `-i PATH` | Input data file |
| `-mesh PATH` | Mesh file, required for MFEM |
| `-fields [F1,F2,...]` | Fields to process. All fields are processed by default |
| `-elements [begin,end)` | Contiguous element range to process |
| `-M ABS\|REL EB` | Error-bound mode and value |
| `-method NAME` | `FDWZ`, `MRCZ`, or a baseline with its layout, such as `SZ3-4D` or `ZFP-3D-fused` |
| `-max-cpus N` | Number of OpenMP threads |
| `-partition-num N` | Number of independently compressed partitions |
| `-reconstructed-output PATH` | Writes the reconstructed field during decompression |

### Experiment Driver

`t1_start_c_code.py` reads a JSON configuration, builds the required dependencies, and runs every combination of field, error bound, method, and thread count:

```bash
python3 t1_start_c_code.py configs/all_small_datasets.json
```

A configuration entry looks like this:

```json
{
  "name": "mixlay0",
  "data_path": "third_party/pySEMTools/examples/data/mixlay0.f00001",
  "format": "nekrs",
  "field_list": ["UX", "UY", "P", "T", "S00", "S01"],
  "eb_mode": "REL",
  "eb_list": [0.01, 0.001, 0.0001, 0.00001, 0.000001],
  "method_list": ["MRCZ", "FDWZ", "SZ3-4D"],
  "max_cpus_list": [1]
}
```

### Using the C++ Library

External CMake projects can link the header-only core:

```cmake
find_package(MRCZ CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE mrcz::mrcz)
```

```cpp
#include <mrcz_core/mrcz.hpp>
```

## How It Works

### FDWZ

FDWZ follows the prediction--quantization--encoding pipeline of SZ but keeps all numerical dependencies inside each element.

1. The prediction levels, visiting order, reference positions, and interpolation weights are precomputed from the element shape and the GLL coordinates of each dimension. The two endpoints of each dimension form the initial boundary, and each subsequent level selects the unpredicted node closest to the midpoint of an interval in reference space.
2. The 4 or 8 corner points of a 2D or 3D element are quantized with zero as the prediction. The remaining nodes are predicted by cubic, quadratic, or linear Lagrange interpolation from reconstructed references.
3. Each prediction error is quantized immediately and the reconstructed value replaces the original, so the compressor and decompressor follow the same error propagation.
4. Each partition samples elements at a fixed interval, compares six dimension orders, and uses the selected order for the final cubic prediction. The default sampling interval is 128 with at least 32 sampled elements.
5. Quantization symbols are grouped by corner points and prediction levels and encoded with up to five Huffman trees. All Huffman metadata and payloads are compressed together in one Zstandard frame.
6. Symbols are stored as 32-bit integers. `INT32_MIN` marks values that cannot be represented safely, and the corresponding original values are appended after the Huffman stream.
7. Partitions are compressed independently in parallel with OpenMP. Prediction stencils and GLL weights are precomputed, and the final prediction processes four elements per batch.

### MRCZ

MRCZ uses an element-wise transform representation as the prediction model and encodes pointwise residuals separately.

1. The 1D evaluation matrices of the discrete Chebyshev transform (DChT) and their inverses are constructed from the element shape and the GLL nodes of each dimension.
2. Each element is transformed with a tensor-product forward transform. All coefficients are quantized with `coeff_eb=0.05*abs_eb`, and the quantized values replace the original coefficients.
3. The inverse transform of the quantized coefficients serves as the prediction for every node. The pointwise residual between the original and predicted values is quantized under the user-specified `abs_eb`, which enforces the error bound.
4. The coefficient symbols are reordered from element-major to mode-major order so that the same mode from different elements is stored contiguously. The residual symbols remain element-major.
5. The coefficient and residual streams each use one Huffman tree and are compressed together in one Zstandard frame.
6. The mode-major reordering is performed with cache-aware tiles that cover up to 32 elements and 32 modes. Each OpenMP thread reuses partition-local coefficient and residual buffers.

## Datasets and Baseline Layouts

The submodules include three small NekRS datasets that can be tested directly:

- `third_party/pySEMTools/examples/data/mixlay0.f00001`
- `third_party/pySEMTools/examples/data/tc_channel0.f00001`
- `third_party/fld_pyutils/demos/data/test0.f00001`

The MFEM and Nektar++ datasets can be regenerated with:

```bash
python3 generate_mfem_and_nektar_data.py
```

The script generates two MFEM examples and the Nektar++ cylinder example in VTU format. The two large-scale NekRS datasets used in the paper are not stored in this repository, and their locations are specified in the JSON configurations.

Baselines are called through their unmodified official APIs. For a partition of `E` elements with element shape `[N0,N1,N2]` and `N=N0*N1*N2`, the candidate layouts are:

| Layout | API call |
|---|---|
| 1D | One call on an `[E*N]` array |
| 2D | One call on an `[E,N]` array |
| fused 3D | One call on an `[E*N0,N1,N2]` array |
| per-element 3D | One call per `[N0,N1,N2]` element, with the length of each bitstream stored |
| 4D | One call on an `[E,N0,N1,N2]` array |
| 4D+GLL | The 4D nonuniform-coordinate interface of MGARD-X |

The baseline layouts used for the large-scale datasets in the paper are SZ3-4D, SZ2-3D-fused, QoZ-4D, ZFP-3D-fused, SPERR-2D, and MGARD-X-1D. SZ2 is excluded from the multicore throughput experiments because its official API modifies mutable process-wide state and cannot be called concurrently. The memory demand of MGARD-X prevents it from running at the same high core counts as the other methods on the large-scale datasets, so the corresponding figures show its measured low-core results together with an ideal linear-scaling projection.

## Reproducing the Experiments and Figures

JSON configurations are stored in `configs/`, and experiment CSV files and intermediate plotting data are stored in `results/`. The main entry points are:

| File | Purpose |
|---|---|
| `t1_start_c_code.py` | Builds and runs the rate--distortion, throughput, and reconstruction experiments |
| `t2_find_three_irregular_elements.py` | Finds three nonplanar NekRS elements that pairwise share faces and generates the Plotly HTML used for the geometry illustration |
| `t3_test_cubic_face_prediction.py` | Compares cross-element and face-parallel cubic prediction on shared faces and summarizes errors and quantized-symbol entropy |
| `generate_mfem_and_nektar_data.py` | Generates the MFEM and Nektar++ test data |
| `py_funcs/summarize_rd.py` | Summarizes per-field rate--distortion results at matched PSNRs |
| `py_funcs/summarize_end_to_end.py` | Combines codec times with Globus transfer records |
| `py_funcs/plot_style.py` | Defines the shared fonts, method order, colors, markers, and axis styles of the paper figures |
| `py_funcs/rasterize_paper_figures.py` | Converts complex source figures to 300 dpi PNG files at their final physical size |

The main paper figures can be regenerated with:

```bash
python3 py_funcs/plot_fdwz_tree_ablation.py
python3 py_funcs/plot_fig8_transform_comparison.py
python3 py_funcs/plot_transform_mode_statistics.py
python3 py_funcs/plot_mrcz_reordering_ablation.py
python3 py_funcs/plot_small_rd.py
python3 py_funcs/plot_large_rd.py
python3 py_funcs/plot_combined_throughput.py
python3 py_funcs/plot_end_to_end.py
python3 py_funcs/plot_transfer_size_time.py
```

The visualization figure also requires reconstructed fields and ParaView slices:

```bash
reconstruct_job=$(sbatch --parsable run_visualization_cr33_reconstruct.sbatch)
sbatch --dependency=afterok:${reconstruct_job} run_visualization_cr33_render.sbatch
```

The rendering job saves 1200 dpi master panels of the original field and the seven methods in `results/visualization/cr33_panel_masters/` and downsamples them to the 300 dpi paper figure. To adjust only the layout, skip reconstruction and slicing:

```bash
python3 py_funcs/plot_pb0_cr33_visualization.py --compose-only
```

Hand-drawn diagrams use their final PNG files directly, and `cross-element.pdf` can be rasterized to a 300 dpi PNG at the final single-column width. Matplotlib figures and field visualizations are generated directly as 300 dpi PNG files at the paper width with 8 pt fonts:

```bash
python3 py_funcs/rasterize_paper_figures.py
```

The paper source is kept in a separate Overleaf clone at `overleaf/`, which is not tracked by this repository. The paper is the authoritative reference for final performance numbers, ablation results, platform details, and statistical definitions.

## Citation

A paper describing FDWZ and MRCZ is currently under review. Citation information will be added after publication.
