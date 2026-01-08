#ifndef __GET_FULL_FIELD_NAME_CPP__
#define __GET_FULL_FIELD_NAME_CPP__

#include "jzzz_param.hpp"
#include <stdexcept>

namespace jzzz
{
    string param_class::get_full_field_name(string field)
    {
        if (field=="X" or field=="x")
            throw runtime_error("还没写坐标压缩");
        if (field=="Y" or field=="y")
            throw runtime_error("还没写坐标压缩");
        if (field=="Z" or field=="z")
            throw runtime_error("还没写坐标压缩");
        if (field=="U" or field=="u")
            return "velocity";
        if (field=="P" or field=="p")
            return "pressure";
        if (field=="T" or field=="t")
            return "temperature";
        if (field[0]=='S' or field[0]=='s')
            return "s"+field.substr(1);
        throw runtime_error("未知的缩写名称："+field);
    }
}

#endif