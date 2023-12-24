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
#define BOHR2A 0.529177210903
#define A2BOHR 1 / BOHR2A
#define BOLTZMANN 3.166811429e-6
#define AU2FS 0.02418884254

inline int nthread = 1;
