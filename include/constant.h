#pragma once

#ifdef _WIN32
#define posix_memalign(p, a, s) (((*(p)) = _aligned_malloc((s), (a))), *(p) ? 0 : errno)
#endif

#if defined(__linux__)
#define OS "LINUX"
#elif defined(_WIN32)
#define OS "WINDOWS"
#else
#define OS "UNKNOWN"
#endif

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#define CFREQ 5140.486777894163
#define BOHR2A 0.529177249
#define A2BOHR 1.889725989

inline int nthread = 1;
