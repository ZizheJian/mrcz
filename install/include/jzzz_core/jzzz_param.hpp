#ifndef __JZZZ_PARAM_H__
#define __JZZZ_PARAM_H__

#include <string>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <cstring>
#include <type_traits>

using namespace std;

namespace jzzz
{
    class header_class
    {
        public:
            uint8_t word_size;
            string realtype;
            vector<uint32_t> element_shape;
            uint32_t element_size;
            uint64_t file_element_num;
            string variables;
            uint32_t passive_scalar_num;
    };

    class block_class
    {
        public:
            uint64_t block_id;
            uint64_t element_num;
            unordered_map<string,string> field_paths;

            uint64_t begin_index,end_index;
    };

    class param_class
    {
        public:
            string header_path;
            string eb_mode;
            float eb;
            bool all_blocks_to_process=false;
            vector<uint64_t> blocks_to_process;
            bool all_fields_to_process=false;
            vector<string> fields_to_process;

            header_class header;
            unordered_map<string,vector<vector<float>>> data_range;
            uint64_t block_num;
            vector<block_class> blocks;

            uint64_t element_num_to_process; // 实际需要压缩的element数量
            vector<float> data;
            vector<vector<unsigned char>> zstd_bs;
            vector<float> decompressed_data;

            float param_eb_scale=0.05;
            float param_eb;
            float zero_f_qb_bound=3;
            vector<float> x;
            vector<float> L,neg_L;
            vector<float> element_base_weight;

            vector<float> simple_invert_matrix(vector<float>& A,uint32_t n);
            vector<float> simple_matrix_multiply(vector<float>& A,vector<float>& B,uint32_t n);
            void read_header_file();
            void read_data_range();
            void read_block_info();
            string get_full_field_name(string field);
            float get_data_range(string field);
            void compress();
            void decompress();
            void quantize_and_overwrite_number(float &x,int64_t &y,float eb,float neg_2_eb);
            void quantize_err(float *pred,float *truth,int64_t *qb,uint32_t size,float eb,float neg_2_eb);
            void quantize_and_overwrite2(float *data,int64_t *qb,uint32_t size,float eb,float neg_2_eb,float zero_freq_bound);

            param_class():
                x{-1,-0.91953391,-0.73877387,-0.47792495,-0.16527896,0.16527896,0.47792495,0.73877387,0.91953391,1}
            {
                // if ((!std::is_same<float,float>::value) && (!std::is_same<float,double>::value))
                //     throw runtime_error("仅能以float或double类型处理数据");
                uint32_t n=x.size();
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
                element_base_weight=vector<float>(8*n*n*n,0);
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
            }
    };
}

#endif