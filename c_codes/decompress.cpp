#include <cstdio>
#include <exception>

#include "mrcz_decompress.hpp"
#include "mrcz_core/parse_decompress_args.hpp"

int main(int argc,char **argv)
{
    try
    {
        auto params=mrcz::parse_decompress_args(argc,argv);
        mrcz::decompress(params);
        return 0;
    }
    catch (exception &error)
    {
        fprintf(stderr,"Error: %s\n",error.what());
        return 1;
    }
}
