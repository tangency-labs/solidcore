// SPDX-License-Identifier: Apache-2.0
// Shared native-call serialization; never define a per-translation-unit mutex.
#include "sc_shim.h"

std::recursive_mutex& sc_native_mutex()
{
    static std::recursive_mutex mutex;
    return mutex;
}

#if defined(__APPLE__) && !defined(__clang__)
#include <math.h>

// macOS with GNU g++ (harness/cpp/build.sh): make sin/cos fusion irrelevant to results.
//
// Both compilers fuse a sin(x)/cos(x) pair into one call, at sites of their own
// choosing: on Linux both call glibc sincos, whose results equal sin/cos (#357 aside).
// On macOS LLVM (the Rust side) calls Apple's __sincos_stret, which differs from Apple's
// sin/cos for some arguments, and GCC, finding no sincos on Darwin, calls cexp(0 + x*i),
// which differs again. Upstream SMLib never calls either itself. These definitions take
// precedence over libSystem's at link time, so every fused pair on both sides of the
// difftest binary returns exactly what separate sin and cos calls return. The volatile
// pointers keep GCC from fusing the two calls here back into cexp.
static double (*volatile sc_sin)(double) = sin;
static double (*volatile sc_cos)(double) = cos;

extern "C" struct __double2 __sincos_stret(double x)
{
    struct __double2 r;
    r.__sinval = sc_sin(x);
    r.__cosval = sc_cos(x);
    return r;
}

extern "C" double _Complex cexp(double _Complex z)
{
    const double m = exp(__real__ z);
    const double y = __imag__ z;
    double _Complex r;
    __real__ r = m * sc_cos(y);
    __imag__ r = m * sc_sin(y);
    return r;
}
#endif
