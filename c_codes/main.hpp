#pragma once

#include <tuple>
#include <string>
#include "jzzz_core/jzzz_param.hpp"

using namespace std;

template<typename T>
void compress(jzzz::param_class<T> &params,string field);
template<typename T>
tuple<double,double,double,double,double> decompress(jzzz::param_class<T> &params,string field);

#include "compress.tpp"
#include "decompress.tpp"