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
#include "huffman.hpp"
#include "../third_party/zstd/lib/zstd.h"

using namespace std;

namespace mrcz
{
    class header_class
    {
        public:
            string original_header;

            uint8_t word_size;
            vector<uint32_t> element_shape;
            uint32_t element_size;
            uint64_t file_element_num;
            string variables;
            uint32_t passive_scalar_num;

            char endian_mode; // '<'表示小端，'>'表示大端
    };

    class block_class
    {
        public:
            uint64_t block_id;
            uint64_t element_num;
            uint64_t begin_index,end_index;//block包含的element在参与压缩的element中的相对位置

            vector<unsigned char> zstd_bs;
    };

    template<typename T>
    vector<T> compute_GLL_nodes(uint32_t n);

    template<typename T> // 表示压缩时使用的类型，可能和原始数据类型不同
    class param_class
    {
        public:
            bool if_chebyshev=true;
            bool if_corner_based_prediction=false;
            bool if_zero_freq_expanding=true;
            bool if_coefficient_transpose=true;
            bool if_zero_map=true;

            string data_path;
            string eb_mode;
            T input_eb;
            unordered_map<string,T> abs_eb,rel_eb;
            string method;
            bool all_fields_to_process=false;
            vector<string> all_fields_in_file;
            vector<string> fields_to_process;
            unordered_map<string,uint64_t> field_begin_address;//field中参与压缩的element在文件中的起始位置
            unordered_map<string,uint64_t> field_end_address;
            bool got_element_range=false;
            uint64_t element_begin_index=0;
            uint64_t element_end_index=0;

            header_class header;
            unordered_map<string,T> data_range_min;
            unordered_map<string,T> data_range_max;
            vector<block_class> blocks;

            uint64_t element_num_to_process; // 实际需要压缩的element数量
            vector<T> data;
            vector<T> decompressed_data;

            int max_cpus=128;
            T coeff_eb_scale=0.05;
            T zero_coeff_qb_bound=3;
            vector<vector<T>> x;
            vector<vector<T>> L,neg_L;
            vector<T> element_base_weight;
            vector<uint32_t> no_corner_pos;

            vector<T> simple_invert_matrix(vector<T>& A,uint32_t n);
            vector<T> simple_matrix_multiply(vector<T>& A,vector<T>& B,uint32_t n);
            void read_header(string header_path);
            void read_data_range();
            void read_block_info();
            void save_header_and_element_map_compress();
            void save_header_and_element_map_decompress();
            void get_max_cpus();
            void save_tail_metadata_compress();
            void save_tail_metadata_decompress();
            T get_data_range(string field);

            void compress(string field);
            void quantize_and_overwrite_corner_numbers(int64_t *corner_qb,uint64_t begin_index,uint64_t end_index,uint32_t n,uint32_t n2,uint32_t n3,
                T abs_eb,T neg_2_eb);
            void subtract_element_base(T *element,uint32_t n,uint32_t n2,uint32_t n3);
            void forward_transform(vector<T> &temp_element1,vector<T> &temp_element2,T *element,uint32_t n0,uint32_t n1,uint32_t n2,uint32_t n12,uint32_t n012,
                vector<T> &neg_L0,vector<T> &neg_L1,vector<T> &neg_L2);
            void inverse_transform(vector<T> &temp_element1,vector<T> &temp_element2,uint32_t n0,uint32_t n1,uint32_t n2,uint32_t n12,uint32_t n012,
                vector<T> &L0,vector<T> &L1,vector<T> &L2);
            void convert_err_eb_to_no_corner_err_eb(int64_t *element_no_corner_err_qb,int64_t *element_err_qb,uint32_t n,uint32_t n2,uint32_t n3);
            void lossless_compress_coeff_qb(uint64_t begin_index,uint64_t end_index,int64_t *coeff_qb,unsigned char *&huffman_bs,uint32_t n012);
            void lossless_compress_coeff_qb_with_zero_map(uint64_t begin_index,uint64_t end_index,int64_t *f_qb,unsigned char *&huffman_bs,uint32_t n012);

            void decompress(string field);
            void lossless_decompress_coeff_qb_with_zero_map(uint64_t begin_index,uint64_t end_index,int64_t *coeff_qb,const unsigned char *&huffman_bs,uint32_t n012);
            void lossless_decompress_coeff_qb(uint64_t begin_index,uint64_t end_index,int64_t *coeff_qb,const unsigned char *&huffman_bs,uint32_t n012);
            void reconstruct_corner_number(int64_t *corner_qb,uint64_t begin_index,uint64_t end_index,uint32_t n,uint32_t n2,uint32_t n3,T abs_eb);
            void add_element_base(T *element,uint32_t n,uint32_t n2,uint32_t n3);

            void quantize_and_overwrite(T &x,int64_t &y,T eb,T neg_2_eb);
            void quantize_and_overwrite_err(T &pred,T &truth,int64_t &qb,T eb,T neg_2_eb);
            void quantize_err_array(T *pred,T *truth,int64_t *qb,uint32_t size,T eb,T neg_2_eb);
            void quantize_and_overwrite_array1(T *data,int64_t *qb,uint32_t size,T eb,T neg_2_eb);
            void quantize_and_overwrite_array2(T *data,int64_t *qb,uint32_t size,T eb,T neg_2_eb,T zero_coeff_bound);

            void default_init()
            {
                if ((!std::is_same<T,float>::value) && (!std::is_same<T,double>::value))
                    throw runtime_error("仅能以float或double类型处理数据");
                if (header.element_shape.size()!=3)
                    throw runtime_error("element必须是三维的");
                for (int dim=0;dim<3;dim++)
                {
                    x.push_back(compute_GLL_nodes<T>(header.element_shape[dim]));
                    uint32_t n=x[dim].size();
                    L.push_back(vector<T>());
                    for (uint32_t i=0;i<n;i++)
                    {
                        if (i==0)
                            for (uint32_t j=0;j<n;j++)
                                L[dim].push_back(1);
                        else if (i==1)
                            for (uint32_t j=0;j<n;j++)
                                L[dim].push_back(x[dim][j]);
                        else
                            for (uint32_t j=0;j<n;j++)
                                if (if_chebyshev)
                                    L[dim].push_back(2*x[dim][j]*L[dim][(i-1)*n+j]-L[dim][(i-2)*n+j]);
                                else
                                    L[dim].push_back(((2*i-1)*x[dim][j]*L[dim][(i-1)*n+j]-(i-1)*L[dim][(i-2)*n+j])/i);
                    }
                    neg_L.push_back(vector<T>());
                    neg_L[dim]=simple_invert_matrix(L[dim],n);
                }
                uint32_t n0=header.element_shape[0];
                uint32_t n1=header.element_shape[1];
                uint32_t n2=header.element_shape[2];
                uint32_t n12=n1*n2;
                uint32_t n012=n0*n1*n2;
                element_base_weight=vector<T>(8*n012,0);
                for (uint32_t i0=0;i0<n0;i0++)
                    for (uint32_t i1=0;i1<n1;i1++)
                        for (uint32_t i2=0;i2<n2;i2++)
                        {
                            uint32_t bias=8*(i0*n12+i1*n2+i2);
                            element_base_weight[bias+0]=((1-x[0][i0])*(1-x[1][i1])*(1-x[2][i2]))/8;
                            element_base_weight[bias+1]=((1-x[0][i0])*(1-x[1][i1])*(x[2][i2]+1))/8;
                            element_base_weight[bias+2]=((1-x[0][i0])*(x[1][i1]+1)*(1-x[2][i2]))/8;
                            element_base_weight[bias+3]=((1-x[0][i0])*(x[1][i1]+1)*(x[2][i2]+1))/8;
                            element_base_weight[bias+4]=((x[0][i0]+1)*(1-x[1][i1])*(1-x[2][i2]))/8;
                            element_base_weight[bias+5]=((x[0][i0]+1)*(1-x[1][i1])*(x[2][i2]+1))/8;
                            element_base_weight[bias+6]=((x[0][i0]+1)*(x[1][i1]+1)*(1-x[2][i2]))/8;
                            element_base_weight[bias+7]=((x[0][i0]+1)*(x[1][i1]+1)*(x[2][i2]+1))/8;
                        }
                for (uint32_t i=1;i<n2-1;i++)
                    no_corner_pos.push_back(i);
                for (uint32_t i=n2;i<n12-n2;i++)
                    no_corner_pos.push_back(i);
                for (uint32_t i=n12-n2+1;i<n12-1;i++)
                    no_corner_pos.push_back(i);
                for (uint32_t i=n12;i<n012-n12;i++)
                    no_corner_pos.push_back(i);
                for (uint32_t i=n012-n12+1;i<n012-n12+n2-1;i++)
                    no_corner_pos.push_back(i);
                for (uint32_t i=n012-n12+n2;i<n012-n2;i++)
                    no_corner_pos.push_back(i);
                for (uint32_t i=n012-n2+1;i<n012-1;i++)
                    no_corner_pos.push_back(i);
            }

            param_class(){}

            template<typename T2>
            param_class(param_class<T2> &other)
            {
                if_chebyshev=other.if_chebyshev;
                if_corner_based_prediction=other.if_corner_based_prediction;
                if_zero_freq_expanding=other.if_zero_freq_expanding;
                if_coefficient_transpose=other.if_coefficient_transpose;
                if_zero_map=other.if_zero_map;
                
                data_path=other.data_path;
                eb_mode=other.eb_mode;
                input_eb=other.input_eb;
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
                element_begin_index=other.element_begin_index;
                element_end_index=other.element_end_index;

                header=other.header;
                for (auto pair:other.data_range_min)
                    data_range_min[pair.first]=pair.second;
                for (auto pair:other.data_range_max)
                    data_range_max[pair.first]=pair.second;
                blocks=other.blocks;

                element_num_to_process=other.element_num_to_process;
                for (auto val:other.data)
                    data.push_back(val);
                for (auto val:other.decompressed_data)
                    decompressed_data.push_back(val);
                
                max_cpus=other.max_cpus;
                coeff_eb_scale=other.coeff_eb_scale;
                zero_coeff_qb_bound=other.zero_coeff_qb_bound;
                default_init();
            }
    };
}

#include "compress.tpp"
#include "decompress.tpp"
#include "get_max_cpus.tpp"
#include "gll_nodes.tpp"
#include "quantize.tpp"
#include "read_header_file.tpp"
#include "save_header_and_element_map.tpp"
#include "save_tail_metadata.tpp"
#include "simple_matrix_tools.tpp"