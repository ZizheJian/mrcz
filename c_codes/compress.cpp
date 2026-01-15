#ifndef __COMPRESS_CPP__
#define __COMPRESS_CPP__

#include <fstream>
#include "jzzz_core/jzzz_param.hpp"

using namespace std;

void compress(jzzz::param_class &params,string field)
{
    string full_field_name=params.get_full_field_name(field);
    if (full_field_name=="pressure")
    {
        params.data.resize(params.element_num_to_process*params.header.element_size);
        string temp_name=full_field_name+"_path";
        for (auto block_id:params.blocks_to_process)
        {
            auto block=params.blocks[block_id];
            string data_path=block.field_paths[temp_name];
            printf("Reading %s data in block %ld\n",full_field_name.c_str(),block_id);
            ifstream file(data_path,ios::binary);
            file.read(reinterpret_cast<char*>(params.data.data()+block.begin_index*params.header.element_size),block.element_num*params.header.element_size*params.header.word_size);
            file.close();
        }
        params.compress();
        for (auto block_id:params.blocks_to_process)
        {
            auto block=params.blocks[block_id];
            string compressed_path=block.field_paths[temp_name]+".jzzz";
            printf("Writing compressed %s data in block %ld\n",full_field_name.c_str(),block_id);
            ofstream file(compressed_path,ios::binary);
            uint64_t zstd_bs_size=params.zstd_bs[block_id].size();
            file.write(reinterpret_cast<char*>(params.zstd_bs[block_id].data()),zstd_bs_size);
            file.close();
            params.zstd_bs[block_id].clear();
            params.zstd_bs[block_id].shrink_to_fit();
        }
        params.data.clear();
        params.data.shrink_to_fit();
        params.zstd_bs.clear();
        params.zstd_bs.shrink_to_fit();
    }
}

#endif