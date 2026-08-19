#pragma once

#include <sstream>
#include <stdexcept>
#include <string>

#include "mrcz.hpp"

namespace mrcz
{
    inline param_class<float> parse_compress_args(int argc,char **argv)
    {
        param_class<float> params;
        for (int i=1;i<argc;i++)
        {
            if (string(argv[i])=="-format")
            {
                if (i+1>=argc)
                    throw runtime_error("缺少-format参数后的数据格式");
                params.format=argv[++i];
            }
            else if (string(argv[i])=="-i")
            {
                if (i+1>=argc)
                    throw runtime_error("缺少-i参数后的输入文件路径");
                params.data_path=argv[++i];
            }
            else if (string(argv[i])=="-mesh")
            {
                if (i+1>=argc)
                    throw runtime_error("缺少-mesh参数后的网格文件路径");
                params.mesh_path=argv[++i];
            }
            else if (string(argv[i])=="-elements")
            {
                if (i+1>=argc)
                    throw runtime_error("缺少-elements参数后的element范围");
                string elements=argv[++i];
                if (elements.size()<5 or elements.front()!='[' or elements.back()!=']' or elements.find(',')==string::npos)
                    throw runtime_error("-elements必须使用[begin,end]格式");
                elements=elements.substr(1,elements.size()-2);
                size_t comma=elements.find(',');
                params.element_begin_index=stoull(elements.substr(0,comma));
                params.element_end_index=stoull(elements.substr(comma+1));
                params.got_element_range=true;
            }
            else if (string(argv[i])=="-fields")
            {
                if (i+1>=argc)
                    throw runtime_error("缺少-fields参数后的字段列表");
                string fields=argv[++i];
                if (fields.size()<2 or fields.front()!='[' or fields.back()!=']')
                    throw runtime_error("-fields必须使用[field1,field2]格式");
                fields=fields.substr(1,fields.size()-2);
                stringstream stream(fields);
                string field;
                while (getline(stream,field,','))
                    if (!field.empty())
                        params.fields_to_process.push_back(field);
            }
            else if (string(argv[i])=="-M")
            {
                if (i+2>=argc)
                    throw runtime_error("缺少-M参数后的误差控制模式和误差值");
                params.eb_mode=argv[++i];
                params.input_eb=stof(argv[++i]);
            }
            else if (string(argv[i])=="-method")
            {
                if (i+1>=argc)
                    throw runtime_error("缺少-method参数后的压缩方法");
                params.method=argv[++i];
            }
            else if (string(argv[i])=="-max-cpus")
            {
                if (i+1>=argc)
                    throw runtime_error("缺少-max-cpus参数后的CPU数量");
                params.max_cpus=stoi(argv[++i]);
                if (params.max_cpus<=0)
                    throw runtime_error("-max-cpus必须是正整数");
            }
            else
                throw runtime_error("未知的命令行参数: "+string(argv[i]));
        }
        return params;
    }
}
