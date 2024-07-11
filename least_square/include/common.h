#pragma once



#if defined(__linux__)
#include "/home/gxt_kt/Projects/debugstream/debugstream.hpp"
#elif defined(__APPLE__)
#include "/Users/gxt_kt/Projects/debugstream/debugstream.hpp"
#else
    #error "Unsupported operating system"
#endif
