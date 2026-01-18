#pragma once

#include <tuple>
#include <string>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include "jzzz_core/jzzz_param.hpp"
#include "sz3_impl/sz3_2_compress.cpp"
#include "sz3_impl/sz3_2_decompress.cpp"

using namespace std;

template<typename T>
void compress(jzzz::param_class<T> &params,string field);
template<typename T>
tuple<double,double,double,double,double> decompress(jzzz::param_class<T> &params,string field);

#include "compress.tpp"
#include "decompress.tpp"