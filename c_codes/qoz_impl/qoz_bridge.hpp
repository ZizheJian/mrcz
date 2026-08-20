#pragma once

#include <cstddef>
#include <vector>

void qoz_compress_4d(const float *data,size_t n0,size_t n1,size_t n2,size_t element_num,
    double abs_eb,std::vector<unsigned char> &compressed);
void qoz_compress_4d(const double *data,size_t n0,size_t n1,size_t n2,size_t element_num,
    double abs_eb,std::vector<unsigned char> &compressed);

void qoz_decompress_4d(const std::vector<unsigned char> &compressed,float *data);
void qoz_decompress_4d(const std::vector<unsigned char> &compressed,double *data);
