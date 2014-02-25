
#include "memory.h"

void * malloc_simd(const size_t size)
{
#if defined WIN32           // WIN32
    return _aligned_malloc(size, 16);
#elif defined __linux__     // Linux
    return memalign(16, size);
#elif defined __MACH__      // Mac OS X
    return malloc(size);
#else                       // other (use valloc for page-aligned memory)
    return valloc(size);
#endif
}

void free_simd(void *data)
{
#if defined WIN32           // WIN32
    _aligned_free(data);
#elif defined __linux__     // Linux
    free(data);
#elif defined __MACH__      // Mac OS X
    free(data);
#else                       // other (use valloc for page-aligned memory)
    free(data);
#endif
}
