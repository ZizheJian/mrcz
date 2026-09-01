#pragma once

#include "../baseline_bitstream.hpp"
#include "../third_party/zfp/include/zfp.h"

template<typename T>
zfp_field *make_zfp_api_field(mrcz::param_class<T> &params,mrcz::partition_class &partition,string method,T *data)
{
    auto &header=params.header;
    uint64_t n0=header.c_element_shape[0];
    uint64_t n1=header.c_element_shape[1];
    uint64_t n2=header.c_element_shape[2];
    zfp_type type=is_same<T,float>::value?zfp_type_float:zfp_type_double;
    if (method=="ZFP-1D")
        return zfp_field_1d(data,type,partition.element_num*header.element_size);
    if (method=="ZFP-2D")
        return zfp_field_2d(data,type,header.element_size,partition.element_num);
    if (method=="ZFP-3D-fused")
        return zfp_field_3d(data,type,n2,n1,partition.element_num*n0);
    if (method=="ZFP-4D")
        return zfp_field_4d(data,type,n2,n1,n0,partition.element_num);
    if (method=="ZFP-3D-per-element")
        return zfp_field_3d(data,type,n2,n1,n0);
    throw runtime_error("未知的ZFP布局: "+method);
}

template<typename T>
vector<unsigned char> zfp_compress_field(mrcz::param_class<T> &params,mrcz::partition_class &partition,string method,T *data,T abs_eb)
{
    zfp_stream *zfp=zfp_stream_open(nullptr);
    zfp_stream_set_accuracy(zfp,abs_eb);
    zfp_field *zfp_data=make_zfp_api_field(params,partition,method,data);
    size_t max_size=zfp_stream_maximum_size(zfp,zfp_data);
    vector<unsigned char> output(max_size);
    bitstream *output_bs=stream_open(output.data(),max_size);
    zfp_stream_set_bit_stream(zfp,output_bs);
    size_t compressed_size=zfp_compress(zfp,zfp_data);
    if (compressed_size==0)
        throw runtime_error("ZFP压缩失败");
    output.resize(compressed_size);
    stream_close(output_bs);
    zfp_field_free(zfp_data);
    zfp_stream_close(zfp);
    return output;
}

template<typename T>
void zfp_compress(mrcz::param_class<T> &params,string field)
{
    auto &header=params.header;
    auto &partitions=params.partitions;
    auto &data=params.data;
    T abs_eb=params.abs_eb[field];
    double start_time=omp_get_wtime();
    #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
    for (uint64_t partition_id=0;partition_id<partitions.size();partition_id++)
    {
        auto &partition=partitions[partition_id];
        uint64_t begin=partition.begin_index*header.element_size;
        if (params.method=="ZFP-3D-per-element")
        {
            partition.zstd_bs.clear();
            for (uint64_t element_id=0;element_id<partition.element_num;element_id++)
            {
                vector<unsigned char> compressed=zfp_compress_field(params,partition,params.method,data.data()+begin+element_id*header.element_size,abs_eb);
                mrcz::baseline::append_uint64(partition.zstd_bs,compressed.size());
                partition.zstd_bs.insert(partition.zstd_bs.end(),compressed.begin(),compressed.end());
            }
        }
        else
            partition.zstd_bs=zfp_compress_field(params,partition,params.method,data.data()+begin,abs_eb);
    }
    double end_time=omp_get_wtime();
    printf("Compression time= %f seconds\n",end_time-start_time);
    printf("\033[91mField= %s, compression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),params.element_num_to_process*header.element_size*header.word_size/((end_time-start_time)*1024*1024*params.max_cpus));
}
