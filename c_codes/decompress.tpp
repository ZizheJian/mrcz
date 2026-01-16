#pragma once

#include <fstream>
#include <cmath>
#include "jzzz_core/jzzz_param.hpp"

template<typename T>
tuple<double,double,double,double,double> decompress(jzzz::param_class<T> &params,string field)
{
    double br=0,cr=0,max_err=0,mse=0,psnr=0;
    params.zstd_bs.resize(params.block_num);
    string temp_name=field+"_path";
    for (auto block_id:params.blocks_to_process)
    {
        auto block=params.blocks[block_id];
        string compressed_path=block.field_paths[temp_name]+".jzzz";
        printf("Reading compressed %s data in block %ld\r",field.c_str(),block_id);
        ifstream file(compressed_path,ios::binary);
        file.seekg(0,ios::end);
        uint64_t file_size=file.tellg();
        file.seekg(0,ios::beg);
        params.zstd_bs[block_id].resize(file_size);
        br+=file_size;
        file.read(reinterpret_cast<char*>(params.zstd_bs[block_id].data()),file_size);
        file.close();
    }
    printf("\n");
    params.decompress();
    br*=(8.0f/(params.element_num_to_process*params.header.element_size));
    cr=(params.header.word_size*8.0f)/br;

    params.data.resize(params.element_num_to_process*params.header.element_size);
    for (auto block_id:params.blocks_to_process)
    {
        auto block=params.blocks[block_id];
        string data_path=block.field_paths[temp_name];
        ifstream file(data_path,ios::binary);
        if (params.header.word_size==4 && sizeof(T)==8)
        {
            vector<float> temp_data(block.element_num*params.header.element_size);
            file.read(reinterpret_cast<char*>(temp_data.data()),block.element_num*params.header.element_size*params.header.word_size);
            for (uint64_t i=0;i<block.element_num*params.header.element_size;i++)
                params.data[block.begin_index*params.header.element_size+i]=temp_data[i];
        }
        else
            file.read(reinterpret_cast<char*>(params.data.data()+block.begin_index*params.header.element_size),block.element_num*params.header.element_size*params.header.word_size);
        file.close();
    }
    for (uint64_t i=0;i<params.element_num_to_process*params.header.element_size;i++)
    {
        T err=params.decompressed_data[i]-params.data[i];
        mse+=err*err;
        if (fabs(err)>max_err)
            max_err=fabs(err);
        if (fabs(err)>params.abs_eb*1.5)
        {
            printf("i=%ld, decompressed= %f, original= %f, err= %f\n",i,params.decompressed_data[i],params.data[i],err);
            getchar();
        }
    }
    mse/=(params.element_num_to_process*params.header.element_size);
    T data_range=params.get_data_range(field);
    psnr=10*log10(data_range*data_range/mse);
    params.data.clear();
    params.data.shrink_to_fit();
    params.decompressed_data.clear();
    params.decompressed_data.shrink_to_fit();
    for (auto &bs:params.zstd_bs)
    {
        bs.clear();
        bs.shrink_to_fit();
    }
    params.zstd_bs.clear();
    params.zstd_bs.shrink_to_fit();
    return make_tuple(br,cr,max_err,mse,psnr);
}