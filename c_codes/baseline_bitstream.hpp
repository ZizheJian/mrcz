#pragma once

#include "mrcz_core/mrcz.hpp"

namespace mrcz::baseline
{
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
