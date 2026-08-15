#pragma once

#include "mrcz_param.hpp"

namespace mrcz
{
    namespace core
    {
        template<typename T>
        void compress(param_class<T> &params);
    }
}

#include "compress.tpp"
