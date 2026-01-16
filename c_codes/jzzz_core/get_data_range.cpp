#ifndef __GET_DATA_RANGE_CPP__
#define __GET_DATA_RANGE_CPP__

#include "jzzz_param.hpp"

namespace jzzz
{
    float param_class::get_data_range(string field)
    {
        if (field=="velocityx")
            return data_range["U"][0][1]-data_range["U"][0][0];
        if (field=="velocityy")
            return data_range["U"][1][1]-data_range["U"][1][0];
        if (field=="velocityz")
            return data_range["U"][2][1]-data_range["U"][2][0];
        if (field=="pressure")
            return data_range["P"][0][1]-data_range["P"][0][0];
        if (field=="temperature")
            return data_range["T"][0][1]-data_range["T"][0][0];
        if ((field[0]=='S' || field[0]=='s') && (field.length()==3))
        {
            int index=stoi(field.substr(1));
            return data_range["S"][index][1]-data_range["S"][index][0];
        }
        throw runtime_error("无法获取指定field的数据范围");
    }
}

#endif