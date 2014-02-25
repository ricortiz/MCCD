#ifndef MCCD_MEMORY_H
#define MCCD_MEMORY_H

#if defined WIN32           // WIN32
  #include <windows.h>
#else
  #include <malloc.h>
#endif

void * malloc_simd(const size_t size);
void free_simd(void *data);

#endif
