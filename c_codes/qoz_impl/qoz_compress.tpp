#pragma once

#include "../baseline_api_layout.hpp"
#include "qoz_bridge.hpp"

template<typename T>
vector<size_t> make_qoz_api_dims(mrcz::param_class<T> &params,mrcz::partition_class &partition,mrcz::baseline::api_layout layout)
{
    auto &header=params.header;
    size_t n0=header.c_element_shape[0];
    size_t n1=header.c_element_shape[1];
    size_t n2=header.c_element_shape[2];
    vector<size_t> dims;
    if (layout==mrcz::baseline::api_layout::flattened_1d)
        dims={partition.element_num*header.element_size};
    else if (layout==mrcz::baseline::api_layout::elements_2d)
        dims={partition.element_num,header.element_size};
    else if (layout==mrcz::baseline::api_layout::fused_3d)
        dims={partition.element_num*n0,n1,n2};
    else if (layout==mrcz::baseline::api_layout::elements_4d)
        dims={partition.element_num,n0,n1,n2};
    else
        dims={n0,n1,n2};
    dims.erase(remove(dims.begin(),dims.end(),1),dims.end());
    if (dims.empty())
        dims.push_back(1);
    return dims;
}

template<typename T>
void qoz_compress(mrcz::param_class<T> &params,string field)
{
    auto &header=params.header;
    auto &partitions=params.partitions;
    auto &data=params.data;
    T abs_eb=params.abs_eb[field];
    double api_abs_eb=abs_eb*0.999;
    auto layout=mrcz::baseline::parse_api_layout(params.method,"QoZ");
    double start_time=omp_get_wtime();
    #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
    for (uint64_t partition_id=0;partition_id<partitions.size();partition_id++)
    {
        auto &partition=partitions[partition_id];
        uint64_t begin=partition.begin_index*header.element_size;
        if (layout==mrcz::baseline::api_layout::per_element_3d)
        {
            partition.zstd_bs.clear();
            mrcz::baseline::append_uint64(partition.zstd_bs,partition.element_num);
            vector<size_t> dims=make_qoz_api_dims(params,partition,layout);
            for (uint64_t element_id=0;element_id<partition.element_num;element_id++)
            {
                vector<unsigned char> compressed;
                qoz_compress_nd(data.data()+begin+element_id*header.element_size,dims,api_abs_eb,compressed);
                mrcz::baseline::append_uint64(partition.zstd_bs,compressed.size());
                partition.zstd_bs.insert(partition.zstd_bs.end(),compressed.begin(),compressed.end());
            }
        }
        else
        {
            vector<size_t> dims=make_qoz_api_dims(params,partition,layout);
            qoz_compress_nd(data.data()+begin,dims,api_abs_eb,partition.zstd_bs);
        }
    }
    double end_time=omp_get_wtime();
    printf("Compression time= %f seconds\n",end_time-start_time);
    printf("\033[91mField= %s, compression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),params.element_num_to_process*header.element_size*header.word_size/((end_time-start_time)*1024*1024*params.max_cpus));
}
