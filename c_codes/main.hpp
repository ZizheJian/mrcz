#pragma once

using namespace std;

#include <tuple>
#include <string>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include "mrcz_core/mrcz_param.hpp"
#include "sz2_impl/sz2_1_compress.cpp"
#include "sz2_impl/sz2_1_decompress.cpp"
#include "sz3_impl/sz3_2_compress.cpp"
#include "sz3_impl/sz3_2_decompress.cpp"
#include "zfp_impl/zfp_compress.cpp"
#include "zfp_impl/zfp_decompress.cpp"
#include "sperr_impl/sperr_compress.cpp"
#include "sperr_impl/sperr_decompress.cpp"
#include "../third_party/get_max_cpus/c_codes/get_max_cpus.hpp"
#include "get_coordinate_in_element.tpp"
#include "mrcz_core/read_and_write_data.tpp"
