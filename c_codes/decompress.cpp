#ifndef __DECOMPRESS_CPP__
#define __DECOMPRESS_CPP__

#include <fstream>
#include <cmath>
#include "jzzz_core/jzzz_param.hpp"

tuple<float,float,float,float> decompress(jzzz::param_class &params,string field)
{
    float br=0,cr=0,mse=0,psnr=0,max_err=0;
    string full_field_name=params.get_full_field_name(field);
    if (full_field_name=="pressure")
    {
        params.zstd_bs.resize(params.block_num);
        string temp_name=full_field_name+"_path";
        for (auto block_id:params.blocks_to_process)
        {
            auto block=params.blocks[block_id];
            string compressed_path=block.field_paths[temp_name]+".jzzz";
            printf("Reading compressed %s data in block %ld\n",full_field_name.c_str(),block_id);
            ifstream file(compressed_path,ios::binary);
            file.seekg(0,ios::end);
            uint64_t file_size=file.tellg();
            file.seekg(0,ios::beg);
            params.zstd_bs[block_id].resize(file_size);
            br+=file_size;
            file.read(reinterpret_cast<char*>(params.zstd_bs[block_id].data()),file_size);
            file.close();
        }
        params.decompress();
        br*=(8.0f/(params.element_num_to_process*params.header.element_size));
        cr=(params.header.word_size*8.0f)/br;

        params.data.resize(params.element_num_to_process*params.header.element_size);
        for (auto block_id:params.blocks_to_process)
        {
            auto block=params.blocks[block_id];
            string data_path=block.field_paths[temp_name];
            ifstream file(data_path,ios::binary);
            file.read(reinterpret_cast<char*>(params.data.data()+block.begin_index*params.header.element_size),block.element_num*params.header.element_size*params.header.word_size);
            file.close();
        }
        for (uint64_t i=0;i<params.element_num_to_process*params.header.element_size;i++)
        {
            float err=params.decompressed_data[i]-params.data[i];
            mse+=err*err;
            if (fabs(err)>max_err)
                max_err=fabs(err);
        }
        mse/=(params.element_num_to_process*params.header.element_size);
        psnr=10*log10((params.data_range[field][0][1]-params.data_range[field][0][0])*(params.data_range[field][0][1]-params.data_range[field][0][0])/mse);
        params.data.clear();
        params.data.shrink_to_fit();
        params.decompressed_data.clear();
        params.decompressed_data.shrink_to_fit();
    }
    printf("max_err= %f\n",max_err);
    return make_tuple(br,cr,mse,psnr);
}

#endif