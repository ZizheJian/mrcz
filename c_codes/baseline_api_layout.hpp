#pragma once

#include "mrcz_core/mrcz.hpp"

namespace mrcz::baseline
{
    enum class api_layout
    {
        flattened_1d,
        elements_2d,
        fused_3d,
        per_element_3d,
        elements_4d
    };

    inline api_layout parse_api_layout(string method,string compressor)
    {
        if (method==compressor+"-1D-uniform")
            return api_layout::flattened_1d;
        if (method==compressor+"-2D-uniform")
            return api_layout::elements_2d;
        if (method==compressor+"-3D-uniform-fused" or method==compressor+"-3D-uniform")
            return api_layout::fused_3d;
        if (method==compressor+"-3D-uniform-per-element")
            return api_layout::per_element_3d;
        if (method==compressor+"-4D-uniform")
            return api_layout::elements_4d;
        throw runtime_error("未知的"+compressor+" API布局: "+method);
    }

    inline void append_uint64(vector<unsigned char> &output,uint64_t value)
    {
        size_t old_size=output.size();
        output.resize(old_size+sizeof(value));
        memcpy(output.data()+old_size,&value,sizeof(value));
    }

    inline uint64_t read_uint64(const vector<unsigned char> &input,size_t &offset)
    {
        if (offset+sizeof(uint64_t)>input.size())
            throw runtime_error("逐element API bitstream的长度表不完整");
        uint64_t value;
        memcpy(&value,input.data()+offset,sizeof(value));
        offset+=sizeof(value);
        return value;
    }
}
