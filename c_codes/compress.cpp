#include <cstdio>
#include <exception>

#include "mrcz_compress.hpp"
#include "mrcz_core/parse_compress_args.hpp"

int main(int argc,char **argv)
{
    try
    {
        auto params=mrcz::parse_compress_args(argc,argv);
        mrcz::compress(params);
        return 0;
    }
    catch (exception &error)
    {
        fprintf(stderr,"Error: %s\n",error.what());
        return 1;
    }
}
