#pragma once

#include "mrcz_param.hpp"

namespace mrcz
{
    namespace core
    {
        template<typename T>
        void decompress(param_class<T> &params);
    }
}

#include "decompress.tpp"
