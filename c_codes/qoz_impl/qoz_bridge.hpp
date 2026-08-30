#pragma once

#include <cstddef>
#include <vector>

void qoz_compress_nd(const float *data,const std::vector<size_t> &dims,double abs_eb,std::vector<unsigned char> &compressed);
void qoz_compress_nd(const double *data,const std::vector<size_t> &dims,double abs_eb,std::vector<unsigned char> &compressed);

void qoz_decompress_nd(const unsigned char *compressed,size_t compressed_size,float *data);
void qoz_decompress_nd(const unsigned char *compressed,size_t compressed_size,double *data);
