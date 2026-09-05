#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <cstring>
#include <type_traits>
#include <omp.h>
#include <immintrin.h>
#include <cmath>
#include <charconv>
#include <limits>
#include <algorithm>
#include <bit>
#include "huffman.hpp"
#include <zstd.h>
#include <mfem.hpp>
#include <tinyxml2.h>

using namespace std;

namespace mrcz
{
    class header_class
    {
        public:
            string original_header;

            uint8_t word_size;
            uint8_t nekrs_space_dimension=3;
            vector<uint32_t> file_element_shape;//文件或原库使用的维度顺序；NekRS为[lx,ly,lz]
            vector<uint32_t> c_element_shape;//当前线性数据按C row-major解释的shape；NekRS为[lz,ly,lx]
            uint32_t element_size;
            uint64_t file_element_num;
            string variables;
            uint32_t passive_scalar_num;

            char endian_mode; // '<'表示小端，'>'表示大端

            uint64_t dof_num=0; //mfem格式中gf文件内dof的总数
            tinyxml2::XMLDocument *header_doc;//用于存储nektar++生成的vtu文件
    };

    class partition_class
    {
        public:
            uint64_t partition_id;
            uint64_t element_num;
            uint64_t begin_index,end_index;//partition包含的element在参与压缩的element中的相对位置

            vector<unsigned char> zstd_bs;
    };

    template<typename T>
    vector<T> compute_GLL_nodes(uint32_t n);

    template<typename T> // 表示压缩时使用的类型，可能和原始数据类型不同
    class param_class
    {
        public:
            string grid_assumption="GLL"; // GLL, chebyshev, uniform
            bool if_constant_compression=true;
            bool if_measure_phase_time=false;

            // MRCZ优化开关
            string transform="DChT"; // DLT, DChT, DST, DCT
            T coeff_eb_scale=0.05;
            bool if_coefficient_transpose=true;
            bool if_tiled_coefficient_transpose=true;

            // FDWZ优化开关
            bool if_max_huffman_tree=true; //仅FDWZ-ablation读取；不写入metadata，解压时必须与压缩一致
            string fdwz_plan="sample"; //仅FDWZ-ablation读取；sample或XYZ/XZY/YXZ/YZX/ZXY/ZYX-linear/cubic
            uint64_t fdwz_sample_interval=128; //0使用旧自动样本数；正整数表示每隔多少个element采一个区间中点

            string format;
            string data_path;
            string mesh_path;
            string eb_mode;
            T input_eb;
            unordered_map<string,T> abs_eb,rel_eb; //每个field的绝对误差和相对误差
            string compute_precision="auto"; //auto、float或double；double输入不会降为float
            string method;
            bool all_fields_to_process=false;
            vector<string> all_fields_in_file;
            vector<string> fields_to_process;
            unordered_map<string,uint64_t> field_begin_address;//field中参与压缩的element在文件中的起始位置
            unordered_map<string,uint64_t> field_end_address;
            bool got_element_range=false;
            uint64_t element_begin_index=0;
            uint64_t element_end_index=0; //内部统一使用[begin,end)半开区间

            header_class header;
            unordered_map<string,T> data_range_min;
            unordered_map<string,T> data_range_max;
            vector<partition_class> partitions;

            uint64_t element_num_to_process; // 实际需要压缩的element数量
            vector<T> data;
            vector<T> decompressed_data;

            int max_cpus=128;
            uint64_t partition_num=0; //0表示压缩时沿用max_cpus；解压时从bitstream读取

            vector<vector<T>> x;
            vector<vector<T>> L,neg_L;

            vector<T> simple_invert_matrix(vector<T>& A,uint32_t n);
            void read_header(string task);
            void read_header_nekrs(string task);
            void read_header_mfem(string task);
            void read_header_nektarpp(string task);
            void save_header_and_element_map(string task);
            void save_header_and_element_map_nekrs(string task);
            void save_header_and_element_map_mfem(string task);
            void save_header_and_element_map_nektarpp(string task);
            void save_tail_metadata_compress();
            void save_tail_metadata_decompress();

            // 压缩器入口
            template<bool measure_phase_time,bool ablation>
            void compress_mrcz(string field);
            template<bool measure_phase_time,bool ablation>
            void decompress_mrcz(string field);
            template<bool measure_phase_time,bool ablation>
            void compress_fdwz(string field);
            template<bool measure_phase_time,bool ablation>
            void decompress_fdwz(string field);

            // MRCZ变换
            void forward_transform(vector<T> &temp_element1,vector<T> &temp_element2,T *element,uint32_t n0,uint32_t n1,uint32_t n2,uint32_t n12,uint32_t n012,
                vector<T> &L0,vector<T> &L1,vector<T> &L2);
            void inverse_transform(vector<T> &temp_element1,vector<T> &temp_element2,uint32_t n0,uint32_t n1,uint32_t n2,uint32_t n12,uint32_t n012,
                vector<T> &neg_L0,vector<T> &neg_L1,vector<T> &neg_L2);

            // MRCZ coefficient QB无损编解码
            template<bool ablation>
            void lossless_compress_coeff_qb(uint64_t begin_index,uint64_t end_index,int64_t *coeff_qb,unsigned char *&huffman_bs,uint32_t n012);
            template<bool ablation>
            void lossless_decompress_coeff_qb(uint64_t begin_index,uint64_t end_index,int64_t *coeff_qb,const unsigned char *&huffman_bs,uint32_t n012);

            // 通用量化
            void quantize_and_overwrite(T &x,int64_t &y,T eb,T neg_2_eb);
            void quantize_and_overwrite_array(T *data,int64_t *qb,uint32_t size,T eb,T neg_2_eb);
            void quantize_and_overwrite_err(T &pred,T &truth,int64_t &qb,T eb,T neg_2_eb);
            void quantize_array(T *pred,T *truth,int64_t *qb,uint32_t size,T eb,T neg_2_eb);

            void default_init()
            {
                if ((!is_same<T,float>::value) and (!is_same<T,double>::value))
                    throw runtime_error("仅能以float或double类型处理数据");
                if (header.c_element_shape.size()!=3)
                    throw runtime_error("element必须是三维的。如果是二维的，应该补上一个长度为1的维度。");
                if (grid_assumption!="GLL" and grid_assumption!="chebyshev" and grid_assumption!="uniform")
                    throw runtime_error("不支持的grid assumption: "+grid_assumption);
                if (transform!="DLT" and transform!="DChT" and transform!="DST" and transform!="DCT")
                    throw runtime_error("不支持的transform: "+transform);
                if (transform=="DST" and grid_assumption!="uniform")
                    throw runtime_error("DST在"+grid_assumption+"网格上是奇异的");

                x.assign(3,{});
                L.assign(3,{});
                neg_L.assign(3,{});
                for (int dim=0;dim<3;dim++)
                {
                    uint32_t n=header.c_element_shape[dim];
                    if (n==0)
                        throw runtime_error("element的每个维度都必须至少包含1个点");
                    if (n==1)
                        x[dim]={0};
                    else if (grid_assumption=="GLL")
                        x[dim]=compute_GLL_nodes<T>(n);
                    else if (grid_assumption=="chebyshev")
                        for (uint32_t i=0;i<n;i++)
                            x[dim].push_back(-cos(M_PI*i/(n-1)));
                    else if (grid_assumption=="uniform")
                        for (uint32_t i=0;i<n;i++)
                            x[dim].push_back(2.0*i/(n-1)-1);

                    T pi=acos(T(-1));
                    for (uint32_t i=0;i<n;i++)
                    {
                        if (transform=="DLT")
                        {
                            if (i==0)
                                for (uint32_t j=0;j<n;j++)
                                    neg_L[dim].push_back(1);
                            else if (i==1)
                                for (uint32_t j=0;j<n;j++)
                                    neg_L[dim].push_back(x[dim][j]);
                            else
                                for (uint32_t j=0;j<n;j++)
                                    neg_L[dim].push_back(((2*i-1)*x[dim][j]*neg_L[dim][(i-1)*n+j]-(i-1)*neg_L[dim][(i-2)*n+j])/i);
                        }
                        else if (transform=="DChT")
                        {
                            if (i==0)
                                for (uint32_t j=0;j<n;j++)
                                    neg_L[dim].push_back(1);
                            else if (i==1)
                                for (uint32_t j=0;j<n;j++)
                                    neg_L[dim].push_back(x[dim][j]);
                            else
                                for (uint32_t j=0;j<n;j++)
                                    neg_L[dim].push_back(2*x[dim][j]*neg_L[dim][(i-1)*n+j]-neg_L[dim][(i-2)*n+j]);
                        }
                        else if (transform=="DST")
                            for(uint32_t j=0;j<n;j++)
                                neg_L[dim].push_back((T(2)/(n+1))*sin(pi*(i+1)*(j+1)/(n+1)));
                        else if (transform=="DCT")
                            for(uint32_t j=0;j<n;j++)
                                neg_L[dim].push_back(cos(pi*i*(x[dim][j]+1)/2));
                    }
                    L[dim]=simple_invert_matrix(neg_L[dim],n);
                }
            }

            param_class(){}

            template<typename T2>
            param_class(param_class<T2> &other,string data_to_copy="")
            {
                if_constant_compression=other.if_constant_compression;
                if_measure_phase_time=other.if_measure_phase_time;
                if_coefficient_transpose=other.if_coefficient_transpose;
                if_tiled_coefficient_transpose=other.if_tiled_coefficient_transpose;
                if_max_huffman_tree=other.if_max_huffman_tree;
                fdwz_plan=other.fdwz_plan;
                fdwz_sample_interval=other.fdwz_sample_interval;
                grid_assumption=other.grid_assumption;
                transform=other.transform;
                
                format=other.format;
                data_path=other.data_path;
                mesh_path=other.mesh_path;
                eb_mode=other.eb_mode;
                input_eb=other.input_eb;
                compute_precision=other.compute_precision;
                for (auto pair:other.abs_eb)
                    abs_eb[pair.first]=pair.second;
                for (auto pair:other.rel_eb)
                    rel_eb[pair.first]=pair.second;
                method=other.method;
                all_fields_to_process=other.all_fields_to_process;
                all_fields_in_file=other.all_fields_in_file;
                fields_to_process=other.fields_to_process;
                field_begin_address=other.field_begin_address;
                field_end_address=other.field_end_address;
                got_element_range=other.got_element_range;
                element_begin_index=other.element_begin_index;
                element_end_index=other.element_end_index;

                header=other.header;
                for (auto pair:other.data_range_min)
                    data_range_min[pair.first]=pair.second;
                for (auto pair:other.data_range_max)
                    data_range_max[pair.first]=pair.second;
                partitions=other.partitions;

                element_num_to_process=other.element_num_to_process;
                if (data_to_copy=="data")
                    data.assign(other.data.begin(),other.data.end());
                else if (data_to_copy=="decompressed_data")
                    decompressed_data.assign(other.decompressed_data.begin(),other.decompressed_data.end());
                else if (!data_to_copy.empty())
                    throw runtime_error("不支持的data_to_copy: "+data_to_copy);

                max_cpus=other.max_cpus;
                partition_num=other.partition_num;
                coeff_eb_scale=other.coeff_eb_scale;
                default_init();
            }
    };

    namespace core
    {
        template<typename T>
        void compress(param_class<T> &params);

        template<typename T>
        void decompress(param_class<T> &params);
    }
}

#include "coefficient_lossless_compress.tpp"
#include "coefficient_lossless_decompress.tpp"
#include "mrcz_compress.tpp"
#include "mrcz_decompress.tpp"
#include "gll_nodes.tpp"
#include "quantize.tpp"
#include "read_header.tpp"
#include "save_header_and_element_map.tpp"
#include "save_tail_metadata.tpp"
#include "simple_matrix_tools.tpp"
#include "compress.tpp"
#include "decompress.tpp"
