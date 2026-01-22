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

namespace jzzz
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
            uint64_t begin_index,end_index;

            vector<unsigned char> zstd_bs;
    };

    template<typename T> // 表示压缩时使用的类型，可能和原始数据类型不同
    class param_class
    {
        public:
            string data_path;
            string eb_mode;
            T input_eb;
            unordered_map<string,T> abs_eb,rel_eb;
            string method;
            bool all_fields_to_process=false;
            vector<string> all_fields_in_file;
            vector<string> fields_to_process;
            unordered_map<string,uint64_t> field_begin_address;
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
            T zero_f_qb_bound=3;
            vector<T> x;
            vector<T> L,neg_L;
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
            T get_data_range(string field);
            void compress(string field);
            void decompress(string field);
            void quantize_and_overwrite(T &x,int64_t &y,T eb,T neg_2_eb);
            void quantize_and_overwrite_err(T &pred,T &truth,int64_t &qb,T eb,T neg_2_eb);
            void quantize_err_array(T *pred,T *truth,int64_t *qb,uint32_t size,T eb,T neg_2_eb);
            void quantize_and_overwrite_array2(T *data,int64_t *qb,uint32_t size,T eb,T neg_2_eb,T zero_freq_bound);

            void default_init()
            {
                if ((!std::is_same<T,float>::value) && (!std::is_same<T,double>::value))
                    throw runtime_error("仅能以float或double类型处理数据");
                uint32_t n=x.size();
                uint32_t n2=n*n;
                uint32_t n3=n2*n;
                for (uint32_t i=0;i<n;i++)
                {
                    if (i==0)
                        for (uint32_t j=0;j<n;j++)
                            L.push_back(1);
                    else if (i==1)
                        for (uint32_t j=0;j<n;j++)
                            L.push_back(x[j]);
                    else
                        for (uint32_t j=0;j<n;j++)
                            L.push_back(2*x[j]*L[(i-1)*n+j]-L[(i-2)*n+j]);
                }
                neg_L=simple_invert_matrix(L,n);
                element_base_weight=vector<T>(8*n*n*n,0);
                for (uint32_t i0=0;i0<n;i0++)
                    for (uint32_t i1=0;i1<n;i1++)
                        for (uint32_t i2=0;i2<n;i2++)
                        {
                            uint32_t bias=8*(i0*n*n+i1*n+i2);
                            element_base_weight[bias+0]=((1-x[i0])*(1-x[i1])*(1-x[i2]))/8;
                            element_base_weight[bias+1]=((1-x[i0])*(1-x[i1])*(x[i2]+1))/8;
                            element_base_weight[bias+2]=((1-x[i0])*(x[i1]+1)*(1-x[i2]))/8;
                            element_base_weight[bias+3]=((1-x[i0])*(x[i1]+1)*(x[i2]+1))/8;
                            element_base_weight[bias+4]=((x[i0]+1)*(1-x[i1])*(1-x[i2]))/8;
                            element_base_weight[bias+5]=((x[i0]+1)*(1-x[i1])*(x[i2]+1))/8;
                            element_base_weight[bias+6]=((x[i0]+1)*(x[i1]+1)*(1-x[i2]))/8;
                            element_base_weight[bias+7]=((x[i0]+1)*(x[i1]+1)*(x[i2]+1))/8;
                        }
                no_corner_pos.resize(n3-8);
                for (uint32_t i=1;i<n-1;i++)
                    no_corner_pos[i-1]=i;
                for (uint32_t i=n;i<n2-n;i++)
                    no_corner_pos[i-2]=i;
                for (uint32_t i=n2-n+1;i<n2-1;i++)
                    no_corner_pos[i-3]=i;
                for (uint32_t i=n2;i<n3-n2;i++)
                    no_corner_pos[i-4]=i;
                for (uint32_t i=n3-n2+1;i<n3-n2+n-1;i++)
                    no_corner_pos[i-5]=i;
                for (uint32_t i=n3-n2+n;i<n3-n;i++)
                    no_corner_pos[i-6]=i;
                for (uint32_t i=n3-n+1;i<n3-1;i++)
                    no_corner_pos[i-7]=i;
            }

            param_class():
                x{-1,-0.91953391,-0.73877387,-0.47792495,-0.16527896,0.16527896,0.47792495,0.73877387,0.91953391,1}
            {
                default_init();
            }

            template<typename T2>
            param_class(param_class<T2> &other)
            {
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
                zero_f_qb_bound=other.zero_f_qb_bound;
                for (auto val:other.x)
                    x.push_back(val);
                default_init();
            }
    };
}

#include "compress.tpp"
#include "decompress.tpp"
#include "get_max_cpus.tpp"
#include "quantize.tpp"
#include "read_header_file.tpp"
#include "save_header_and_element_map.tpp"
#include "simple_matrix_tools.tpp"