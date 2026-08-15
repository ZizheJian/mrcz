#pragma once

#include "parse_compress_args.hpp"

namespace mrcz
{
    inline param_class<float> parse_decompress_args(int argc,char **argv)
    {
        return parse_compress_args(argc,argv);
    }
}
