#ifndef COMMON_HH
#define COMMON_HH

#include <cstdio>
#include <cstdlib>

FILE* safe_fopen(const char* filename,const char* mode);

void fatal_error(const char *p,int code);

#endif
