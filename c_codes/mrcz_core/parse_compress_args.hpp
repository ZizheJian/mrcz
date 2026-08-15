#pragma once

#include <sstream>
#include <stdexcept>
#include <string>

#include "mrcz_param.hpp"

namespace mrcz
{
    inline param_class<float> parse_compress_args(int argc,char **argv)
    {
        param_class<float> params;
        for (int i=1;i<argc;i++)
        {
            const std::string arg=argv[i];
            if (arg=="-format")
            {
                if (i+1>=argc)
                    throw std::runtime_error("缺少-format参数后的数据格式");
                params.format=argv[++i];
            }
            else if (arg=="-i")
            {
                if (i+1>=argc)
                    throw std::runtime_error("缺少-i参数后的输入文件路径");
                params.data_path=argv[++i];
            }
            else if (arg=="-mesh")
            {
                if (i+1>=argc)
                    throw std::runtime_error("缺少-mesh参数后的网格文件路径");
                params.mesh_path=argv[++i];
            }
            else if (arg=="-elements")
            {
                if (i+1>=argc)
                    throw std::runtime_error("缺少-elements参数后的element范围");
                std::string elements=argv[++i];
                if (elements.size()<5 or elements.front()!='[' or elements.back()!=']' or elements.find(',')==std::string::npos)
                    throw std::runtime_error("-elements必须使用[begin,end]格式");
                elements=elements.substr(1,elements.size()-2);
                const size_t comma=elements.find(',');
                params.element_begin_index=std::stoull(elements.substr(0,comma));
                params.element_end_index=std::stoull(elements.substr(comma+1));
                params.got_element_range=true;
            }
            else if (arg=="-fields")
            {
                if (i+1>=argc)
                    throw std::runtime_error("缺少-fields参数后的字段列表");
                std::string fields=argv[++i];
                if (fields.size()<2 or fields.front()!='[' or fields.back()!=']')
                    throw std::runtime_error("-fields必须使用[field1,field2]格式");
                fields=fields.substr(1,fields.size()-2);
                std::stringstream stream(fields);
                std::string field;
                while (std::getline(stream,field,','))
                    if (!field.empty())
                        params.fields_to_process.push_back(field);
            }
            else if (arg=="-M")
            {
                if (i+2>=argc)
                    throw std::runtime_error("缺少-M参数后的误差控制模式和误差值");
                params.eb_mode=argv[++i];
                params.input_eb=std::stof(argv[++i]);
            }
            else if (arg=="-method")
            {
                if (i+1>=argc)
                    throw std::runtime_error("缺少-method参数后的压缩方法");
                params.method=argv[++i];
            }
            else if (arg=="-max-cpus")
            {
                if (i+1>=argc)
                    throw std::runtime_error("缺少-max-cpus参数后的CPU数量");
                params.max_cpus=std::stoi(argv[++i]);
                if (params.max_cpus<=0)
                    throw std::runtime_error("-max-cpus必须是正整数");
            }
            else
                throw std::runtime_error("未知的命令行参数: "+arg);
        }
        return params;
    }
}
