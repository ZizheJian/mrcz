#include "qoz_bridge.hpp"

#include "QoZ/api/sz.hpp"

using namespace std;

template<typename T>
void compress_nd(const T *data,const vector<size_t> &dims,double abs_eb,vector<unsigned char> &compressed)
{
    QoZ::Config config;
    config.setDims(dims.begin(),dims.end());
    config.errorBoundMode=QoZ::EB_ABS;
    config.absErrorBound=abs_eb;
    config.QoZ=4;
    config.openmp=false;
    config.verbose=false;
    if (dims.size()==1)
        config.cmprAlgo=QoZ::ALGO_INTERP;

    size_t compressed_size=0;
    char *compressed_data=SZ_compress<T>(config,data,compressed_size);
    compressed.assign(compressed_data,compressed_data+compressed_size);
    delete[] compressed_data;
}

template<typename T>
void decompress_nd(const unsigned char *compressed,size_t compressed_size,T *data)
{
    QoZ::Config config;
    T *output=data;
    SZ_decompress<T>(config,reinterpret_cast<char *>(const_cast<unsigned char *>(compressed)),compressed_size,output);
}

void qoz_compress_nd(const float *data,const vector<size_t> &dims,double abs_eb,vector<unsigned char> &compressed)
{
    compress_nd(data,dims,abs_eb,compressed);
}

void qoz_compress_nd(const double *data,const vector<size_t> &dims,double abs_eb,vector<unsigned char> &compressed)
{
    compress_nd(data,dims,abs_eb,compressed);
}

void qoz_decompress_nd(const unsigned char *compressed,size_t compressed_size,float *data)
{
    decompress_nd(compressed,compressed_size,data);
}

void qoz_decompress_nd(const unsigned char *compressed,size_t compressed_size,double *data)
{
    decompress_nd(compressed,compressed_size,data);
}
