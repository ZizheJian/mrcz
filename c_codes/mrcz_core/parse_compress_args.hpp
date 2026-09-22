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
        if (elements.size()<5 or elements.front()!='[' or (elements.back()!=']' and elements.back()!=')'))
            throw runtime_error("-elements必须使用[begin,end]或[begin,end)格式");

        char right_bracket=elements.back();
        string range_body=elements.substr(1,elements.size()-2);
        size_t comma=range_body.find(',');
        if (comma==string::npos or comma==0 or comma==range_body.size()-1 or range_body.find(',',comma+1)!=string::npos)
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
            else if (string(argv[i])=="-constant-compression")
            {
                if (i+1>=argc)
                    throw runtime_error("缺少-constant-compression参数后的布尔值");
                params.if_constant_compression=parse_bool(argv[++i],"-constant-compression");
            }
            else if (string(argv[i])=="-measure-phase-time")
            {
                if (i+1>=argc)
                    throw runtime_error("缺少-measure-phase-time参数后的布尔值");
                params.if_measure_phase_time=parse_bool(argv[++i],"-measure-phase-time");
            }
            else if (string(argv[i])=="-grid-assumption")
            {
                if (i+1>=argc)
                    throw runtime_error("缺少-grid-assumption参数后的网格假设");
                params.grid_assumption=argv[++i];
                if (params.grid_assumption!="GLL" and params.grid_assumption!="chebyshev" and params.grid_assumption!="uniform")
                    throw runtime_error("-grid-assumption必须是GLL、chebyshev或uniform");
            }
            else if (string(argv[i])=="-transform")
            {
                if (i+1>=argc)
                    throw runtime_error("缺少-transform参数后的变换名称");
                params.transform=argv[++i];
                if (params.transform!="DLT" and params.transform!="DChT" and params.transform!="DST" and params.transform!="DCT")
                    throw runtime_error("-transform必须是DLT、DChT、DST或DCT");
            }
            else if (string(argv[i])=="-coeff-eb-scale")
            {
                if (i+1>=argc)
                    throw runtime_error("缺少-coeff-eb-scale参数后的正数");
                params.coeff_eb_scale=parse_positive_float(argv[++i],"-coeff-eb-scale");
            }
            else if (string(argv[i])=="-coefficient-transpose")
            {
                if (i+1>=argc)
                    throw runtime_error("缺少-coefficient-transpose参数后的布尔值");
                params.if_coefficient_transpose=parse_bool(argv[++i],"-coefficient-transpose");
            }
            else if (string(argv[i])=="-tiled-coefficient-transpose")
            {
                if (i+1>=argc)
                    throw runtime_error("缺少-tiled-coefficient-transpose参数后的布尔值");
                params.if_tiled_coefficient_transpose=parse_bool(argv[++i],"-tiled-coefficient-transpose");
            }
            else if (string(argv[i])=="-max-huffman-tree")
            {
                if (i+1>=argc)
                    throw runtime_error("缺少-max-huffman-tree参数后的布尔值");
                params.if_max_huffman_tree=parse_bool(argv[++i],"-max-huffman-tree");
            }
            else if (string(argv[i])=="-analyze-fdwz-prediction-locality")
            {
                if (i+1>=argc)
                    throw runtime_error("缺少-analyze-fdwz-prediction-locality参数后的布尔值");
                params.if_analyze_fdwz_prediction_locality=parse_bool(argv[++i],"-analyze-fdwz-prediction-locality");
            }
            else if (string(argv[i])=="-compute-precision")
            {
                if (i+1>=argc)
                    throw runtime_error("缺少-compute-precision参数后的计算精度");
                params.compute_precision=argv[++i];
                if (params.compute_precision!="auto" and params.compute_precision!="float" and params.compute_precision!="double")
                    throw runtime_error("-compute-precision必须是auto、float或double");
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
            else if (string(argv[i])=="-partition-num")
            {
                if (i+1>=argc)
                    throw runtime_error("缺少-partition-num参数后的partition数量");
                params.partition_num=parse_uint64(argv[++i],"-partition-num");
                if (params.partition_num==0)
                    throw runtime_error("-partition-num必须是正整数");
            }
            else
                throw runtime_error("未知的命令行参数: "+string(argv[i]));
        }
        return params;
    }
}
