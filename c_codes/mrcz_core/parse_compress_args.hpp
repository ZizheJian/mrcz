#pragma once

#include <sstream>
#include <stdexcept>
#include <string>
#include <limits>

#include <parsers.hpp>
#include "mrcz.hpp"

namespace mrcz
{
    inline void parse_element_range(string elements,uint64_t &begin,uint64_t &end)
    {
        if (elements.size()<5 or elements.front()!='[' or
            (elements.back()!=']' and elements.back()!=')'))
            throw runtime_error("-elements必须使用[begin,end]或[begin,end)格式");

        char right_bracket=elements.back();
        string range_body=elements.substr(1,elements.size()-2);
        size_t comma=range_body.find(',');
        if (comma==string::npos or comma==0 or comma==range_body.size()-1 or
            range_body.find(',',comma+1)!=string::npos)
            throw runtime_error("-elements必须使用[begin,end]或[begin,end)格式");

        string begin_string=range_body.substr(0,comma);
        string end_string=range_body.substr(comma+1);
        begin=parse_uint64(begin_string,"-elements的begin");
        uint64_t input_end=parse_uint64(end_string,"-elements的end");

        if (right_bracket==']')
        {
            if (input_end==numeric_limits<uint64_t>::max())
                throw runtime_error("-elements闭区间的end过大");
            end=input_end+1;
        }
        else
            end=input_end;
    }

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
                parse_element_range(elements,params.element_begin_index,params.element_end_index);
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
                params.input_eb=parse_positive_float(argv[++i],"-M的error bound");
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
                uint32_t max_cpus=parse_positive_uint32(argv[++i],"-max-cpus");
                if (max_cpus>numeric_limits<int>::max())
                    throw runtime_error("-max-cpus超出int范围");
                params.max_cpus=max_cpus;
            }
            else
                throw runtime_error("未知的命令行参数: "+string(argv[i]));
        }
        return params;
    }
}
