#pragma once

#include "qoz_compress.tpp"

template<typename T>
void qoz_decompress(mrcz::param_class<T> &params,string field)
{
    auto &header=params.header;
    auto &partitions=params.partitions;
    auto &decompressed_data=params.decompressed_data;
    auto layout=mrcz::baseline::parse_api_layout(params.method,"QoZ");
    decompressed_data.resize(params.element_num_to_process*header.element_size);
    double start_time=omp_get_wtime();
    #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
    for (uint64_t partition_id=0;partition_id<partitions.size();partition_id++)
    {
        auto &partition=partitions[partition_id];
        uint64_t begin=partition.begin_index*header.element_size;
        if (partition.zstd_bs.empty())
            throw runtime_error(params.method+" partition bitstream为空");
        if (layout==mrcz::baseline::api_layout::per_element_3d)
        {
            size_t offset=0;
            uint64_t element_num=mrcz::baseline::read_uint64(partition.zstd_bs,offset);
            if (element_num!=partition.element_num)
                throw runtime_error(params.method+"逐element bitstream中的element数量不匹配");
            for (uint64_t element_id=0;element_id<element_num;element_id++)
            {
                uint64_t compressed_size=mrcz::baseline::read_uint64(partition.zstd_bs,offset);
                if (compressed_size>partition.zstd_bs.size()-offset)
                    throw runtime_error(params.method+"逐element bitstream不完整");
                qoz_decompress_nd(partition.zstd_bs.data()+offset,compressed_size,decompressed_data.data()+begin+element_id*header.element_size);
                offset+=compressed_size;
            }
            if (offset!=partition.zstd_bs.size())
                throw runtime_error(params.method+"逐element bitstream包含多余数据");
        }
        else
            qoz_decompress_nd(partition.zstd_bs.data(),partition.zstd_bs.size(),decompressed_data.data()+begin);
    }
    double end_time=omp_get_wtime();
    printf("Decompression time= %f seconds\n",end_time-start_time);
    printf("\033[91mField= %s, decompression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),params.element_num_to_process*header.element_size*header.word_size/((end_time-start_time)*1024*1024*params.max_cpus));
}
