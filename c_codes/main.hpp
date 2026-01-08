#ifndef __MAIN_HPP__
#define __MAIN_HPP__

#include <tuple>
#include <string>
#include "jzzz_core/jzzz_param.hpp"

using namespace std;

void compress(jzzz::param_class &params,string field);
tuple<float,float,float,float> decompress(jzzz::param_class &params,string field);

#endif