#pragma once

#include "../baseline_bitstream.hpp"
#pragma push_macro("REL")
#pragma push_macro("ABS")
#undef REL
#undef ABS
#include "mgard/compress_x.hpp"

namespace mrcz::baseline::mgard_detail
{
    constexpr auto absolute_error_bound=mgard_x::error_bound_type::ABS;
}

#pragma pop_macro("ABS")
#pragma pop_macro("REL")

template<typename T>
struct mgard_api_grid
{
    vector<mgard_x::SIZE> shape;
    vector<vector<T>> coords;
    bool force_single_block=false;
};

inline bool mgard_is_per_element(const string &method)
{
    return method=="MGARD-3D-per-element" or method=="MGARD-3D-per-element-GLL";
}

inline bool mgard_is_4d(const string &method)
{
    return method=="MGARD-4D" or method=="MGARD-4D-GLL";
}

template<typename T>
vector<T> make_mgard_axis_coordinates(uint32_t n,const string &grid_assumption)
{
    vector<T> coordinates;
    if(n==1)
        coordinates={T(0)};
    else if(grid_assumption=="GLL")
        coordinates=mrcz::compute_GLL_nodes<T>(n);
    else if(grid_assumption=="chebyshev")
    {
        coordinates.resize(n);
        for(uint32_t i=0;i<n;i++)
            coordinates[i]=-cos(acos(T(-1))*i/(n-1));
    }
    else if(grid_assumption=="uniform")
    {
        coordinates.resize(n);
        for(uint32_t i=0;i<n;i++)
            coordinates[i]=T(i)/(n-1);
    }
    else
        throw runtime_error("MGARD不支持的grid assumption: "+grid_assumption);
    if constexpr(is_same_v<T,double>)
        for(auto &coordinate:coordinates)
            coordinate=static_cast<double>(static_cast<float>(coordinate));
    return coordinates;
}

template<typename T>
mgard_api_grid<T> make_mgard_api_grid(mrcz::param_class<T> &params,mrcz::partition_class &partition,string method)
{
    auto &header=params.header;
    mgard_api_grid<T> grid;
    vector<mgard_x::SIZE> raw_shape;
    if(method=="MGARD-1D")
        raw_shape={partition.element_num*header.element_size};
    else if(method=="MGARD-2D")
        raw_shape={partition.element_num,header.element_size};
    else if(method=="MGARD-3D-fused")
        raw_shape={partition.element_num*header.c_element_shape[0],header.c_element_shape[1],header.c_element_shape[2]};
    else if(mgard_is_per_element(method))
        raw_shape={header.c_element_shape[0],header.c_element_shape[1],header.c_element_shape[2]};
    else if(mgard_is_4d(method))
        raw_shape={partition.element_num,header.c_element_shape[0],header.c_element_shape[1],header.c_element_shape[2]};
    else
        throw runtime_error("未知的MGARD布局: "+method);

    vector<vector<T>> raw_coords;
    string grid_assumption=(method=="MGARD-3D-per-element-GLL" or method=="MGARD-4D-GLL")?"GLL":params.grid_assumption;
    if(grid_assumption!="uniform")
    {
        if(!mgard_is_per_element(method) and !mgard_is_4d(method))
            throw runtime_error(method+"融合了逻辑维度，不能使用显式"+grid_assumption+"坐标");
        if(mgard_is_4d(method))
            raw_coords.push_back(make_mgard_axis_coordinates<T>(partition.element_num,"uniform"));
        for(uint32_t dim=0;dim<3;dim++)
            raw_coords.push_back(make_mgard_axis_coordinates<T>(header.c_element_shape[dim],grid_assumption));
        if(raw_coords.size()!=raw_shape.size())
            throw runtime_error(method+"坐标维数与数组维数不匹配");
    }

    for(size_t dim=0;dim<raw_shape.size();dim++)
    {
        if(raw_shape[dim]==1)
            continue;
        grid.shape.push_back(raw_shape[dim]);
        if(!raw_coords.empty())
            grid.coords.push_back(move(raw_coords[dim]));
    }
    if(grid.shape.empty())
    {
        grid.shape.push_back(1);
        if(!raw_coords.empty())
            grid.coords.push_back({T(0)});
    }
    grid.force_single_block=!grid.coords.empty();
    return grid;
}

inline mgard_x::Config make_mgard_config()
{
    mgard_x::Config config;
    config.dev_type=mgard_x::device_type::SERIAL;
    config.auto_pin_host_buffers=false;
    return config;
}

template<typename T>
void mgard_compress_array(const T *data,const mgard_api_grid<T> &grid,T abs_eb,vector<unsigned char> &output)
{
    void *compressed_data=nullptr;
    size_t compressed_size=0;
    auto config=make_mgard_config();
    if(grid.force_single_block)
    {
        config.domain_decomposition=mgard_x::domain_decomposition_type::Block;
        config.block_size=*max_element(grid.shape.begin(),grid.shape.end());
    }
    auto dtype=is_same_v<T,float>?mgard_x::data_type::Float:mgard_x::data_type::Double;
    mgard_x::compress_status_type status;
    if(grid.coords.empty())
        status=mgard_x::compress(grid.shape.size(),dtype,grid.shape,abs_eb,numeric_limits<double>::infinity(),mrcz::baseline::mgard_detail::absolute_error_bound,data,compressed_data,compressed_size,config,false);
    else
    {
        vector<const mgard_x::Byte*> coord_ptrs;
        for(const auto &coord:grid.coords)
            coord_ptrs.push_back(reinterpret_cast<const mgard_x::Byte*>(coord.data()));
        status=mgard_x::compress(grid.shape.size(),dtype,grid.shape,abs_eb,numeric_limits<double>::infinity(),mrcz::baseline::mgard_detail::absolute_error_bound,data,compressed_data,compressed_size,coord_ptrs,config,false);
    }
    if(status!=mgard_x::compress_status_type::Success)
    {
        free(compressed_data);
        throw runtime_error("MGARD压缩失败，状态码="+to_string(static_cast<int>(status)));
    }
    output.assign(static_cast<unsigned char*>(compressed_data),static_cast<unsigned char*>(compressed_data)+compressed_size);
    free(compressed_data);
}

template<typename T>
void mgard_compress(mrcz::param_class<T> &params,string field)
{
    auto &header=params.header;
    auto &partitions=params.partitions;
    auto &data=params.data;
    T abs_eb=params.abs_eb[field];
    double start_time=omp_get_wtime();
    #pragma omp parallel for schedule(static) num_threads(params.max_cpus)
    for(uint64_t partition_id=0;partition_id<partitions.size();partition_id++)
    {
        auto &partition=partitions[partition_id];
        uint64_t begin=partition.begin_index*header.element_size;
        auto grid=make_mgard_api_grid(params,partition,params.method);
        if(mgard_is_per_element(params.method))
        {
            partition.zstd_bs.clear();
            for(uint64_t element_id=0;element_id<partition.element_num;element_id++)
            {
                vector<unsigned char> compressed;
                mgard_compress_array(data.data()+begin+element_id*header.element_size,grid,abs_eb,compressed);
                mrcz::baseline::append_uint64(partition.zstd_bs,compressed.size());
                partition.zstd_bs.insert(partition.zstd_bs.end(),compressed.begin(),compressed.end());
            }
        }
        else
            mgard_compress_array(data.data()+begin,grid,abs_eb,partition.zstd_bs);
    }
    double end_time=omp_get_wtime();
    printf("Compression time= %f seconds\n",end_time-start_time);
    printf("\033[91mField= %s, compression bandwidth per cpu= %f MB/s\n\033[0m",field.c_str(),params.element_num_to_process*header.element_size*header.word_size/((end_time-start_time)*1024*1024*params.max_cpus));
}
