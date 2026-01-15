#include "main.hpp"
#include <stdexcept>
#include <fstream>
#include <sstream>

int main(int argc,char **argv)
{
    jzzz::param_class params;
    for (int i=1;i<argc;i++)
    {
        string arg=argv[i];
        if (arg=="-i")
        {
            if (i+1<argc)
                params.header_path=argv[++i];
            else
                throw runtime_error("缺少-i参数后的输入文件路径");
            continue;
        }
        if (arg=="-blocks")
        {
            if (i+1<argc)
            {
                string blocks_str=argv[++i];
                blocks_str=blocks_str.substr(1,blocks_str.size()-2);
                stringstream ss(blocks_str);
                string block_id_str;
                while (getline(ss,block_id_str,','))
                    try
                    {
                        uint64_t block_id=stoull(block_id_str);
                        params.blocks_to_process.push_back(block_id);
                    }
                    catch (const std::exception& e)
                    {
                        throw runtime_error("块编号必须为非负整数");
                    }
            }
            else
                throw runtime_error("缺少-blocks参数后的块编号列表");
            continue;
        }
        if (arg=="-M")
        {
            if (i+2<argc)
            {
                params.eb_mode=argv[++i];
                params.eb=stof(argv[++i]);
            }
            else
                throw runtime_error("缺少-M参数后的误差控制模式和误差值");
            continue;
        }
    }
    
    if (params.header_path.empty())
        throw runtime_error("必须指定header文件路径，使用-i参数");
    if (params.eb_mode.empty())
        throw runtime_error("必须指定误差控制模式，使用-M参数");
    if (params.eb<=0.0f)
        throw runtime_error("误差值必须为正数");
    if (params.blocks_to_process.empty())
        params.all_blocks_to_process=true;
    else
        params.all_blocks_to_process=false;
    if (params.fields_to_process.empty())
        params.all_fields_to_process=true;
    else
        params.all_fields_to_process=false;
    
    params.read_header_file();
    params.read_data_range();
    params.read_block_info();
    
    for (string field:params.fields_to_process)
    {
        if (field=="X")
            continue;
        compress(params,field);
        auto result=decompress(params,field);
        float field_br=get<0>(result);
        float field_cr=get<1>(result);
        float field_mse=get<2>(result);
        float field_psnr=get<3>(result);
        printf("Field= %s , %s_EB= %f , BR= %f , CR= %f , MSE= %f , PSNR= %f dB\n",
               field.c_str(),params.eb_mode.c_str(),params.eb,field_br,field_cr,field_mse,field_psnr);
    }
}