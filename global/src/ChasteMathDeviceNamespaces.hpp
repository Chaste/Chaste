/*

Copyright (c) 2005-2024, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef CHASTEMATHDEVICENAMESPACES_HPP_
#define CHASTEMATHDEVICENAMESPACES_HPP_

#include "MathsCustomFunctions.hpp" // Only used this for Signum, candidate for forward delcaration in future
#include "HostDeviceMacros.hpp"
#include <cmath>
#include <type_traits>
#include <algorithm>

/*
    These namespaces allow the printer in chaste-codegen to be device agnostic by
    allowing the generation of function calls such as CHASTE_MATH::Pow(x, y) instead
    of generating pow(x, y).
    
    For cpu calls, we just resolve the call to a std call or one of the methods defined in
    MathsCustomFunctions.hpp.
    
    For device calls, conveniently AMD and NVIDIA have the same API for math calls so we only
    need one generic "Device" namespace that handles both. Instead of calling std on device
    (which is not officially supported), we can call AMD and NVIDIA's APIs. We perform 
    template argument deduction to allow switching between float and double math on device
    which is much more essential for performance on device than cpu. 
*/

namespace ChasteCpu {
    inline constexpr double E = 2.71828182845904523536;
    inline constexpr double Pi = 3.14159265358979323846;
    inline const double NaN = std::numeric_limits<double>::quiet_NaN();
    template<typename T> inline constexpr T Const(double v) { return static_cast<T>(v); }

    // 1-Argument Wrapper
    #define CPU_MATH_1ARG(FuncName, StdFunc) \
    template<typename T> inline T FuncName(T x) { \
        return std::StdFunc(x); \
    }

    // 2-Argument Wrapper
    #define CPU_MATH_2ARG(FuncName, StdFunc) \
    template<typename T1, typename T2> \
    inline auto FuncName(T1 a, T2 b) { \
        return std::StdFunc(a, b); \
    }

    // Reciprocal Wrapper (1.0 / Function)
    #define CPU_MATH_RECIP(FuncName, BaseFunc) \
    template<typename T> inline T FuncName(T x) { \
        return T(1.0) / std::BaseFunc(x); \
    }

    // Inverse Reciprocal Wrapper (Function(1.0 / x))
    #define CPU_MATH_INV_RECIP(FuncName, StdFunc) \
    template<typename T> inline T FuncName(T x) { \
        return std::StdFunc(T(1.0) / x); \
    }

    CPU_MATH_1ARG(Abs,   abs)
    CPU_MATH_1ARG(Acos,  acos)
    CPU_MATH_1ARG(Acosh, acosh)
    CPU_MATH_1ARG(Asin,  asin)
    CPU_MATH_1ARG(Asinh, asinh)
    CPU_MATH_1ARG(Atan,  atan)
    CPU_MATH_1ARG(Atanh, atanh)
    CPU_MATH_1ARG(Ceil,  ceil)
    CPU_MATH_1ARG(Cos,   cos)
    CPU_MATH_1ARG(Cosh,  cosh)
    CPU_MATH_1ARG(Exp,   exp)
    CPU_MATH_1ARG(Expm1, expm1)
    CPU_MATH_1ARG(Floor, floor)
    CPU_MATH_1ARG(Log,   log)
    CPU_MATH_1ARG(Log10, log10)
    CPU_MATH_1ARG(Log1p, log1p)
    CPU_MATH_1ARG(Log2,  log2)
    CPU_MATH_1ARG(Sin,   sin)
    CPU_MATH_1ARG(Sinh,  sinh)
    CPU_MATH_1ARG(Sqrt,  sqrt)
    CPU_MATH_1ARG(Tan,   tan)
    CPU_MATH_1ARG(Tanh,  tanh)

    CPU_MATH_2ARG(Atan2, atan2)
    CPU_MATH_2ARG(Pow,   pow)

    CPU_MATH_RECIP(Sec,  cos)
    CPU_MATH_RECIP(Csc,  sin)
    CPU_MATH_RECIP(Cot,  tan)
    CPU_MATH_RECIP(Sech, cosh)
    CPU_MATH_RECIP(Csch, sinh)
    CPU_MATH_RECIP(Coth, tanh)

    CPU_MATH_INV_RECIP(Asec,  acos)
    CPU_MATH_INV_RECIP(Acsc,  asin)
    CPU_MATH_INV_RECIP(Acot,  atan)
    CPU_MATH_INV_RECIP(Asech, acosh)
    CPU_MATH_INV_RECIP(Acsch, asinh)
    CPU_MATH_INV_RECIP(Acoth, atanh)

    template<typename T> inline T Sign(T x) {
        return Signum(x);
    }

    #undef CPU_MATH_1ARG
    #undef CPU_MATH_2ARG
    #undef CPU_MATH_RECIP
    #undef CPU_MATH_INV_RECIP
}

// Use of this namespace should always be guarded by if USING_DEVICE_COMPILER, as these functions
// are only intended to be inlined into device kernels and other use is untested. 
namespace ChasteDevice {
    // the codegen printer handles wrapping these literals with CHASTE_CONST so they are converted to float if necessary
    inline constexpr double E = 2.71828182845904523536;
    inline constexpr double Pi = 3.14159265358979323846;
#if IN_DEVICE_PASS
    // ChasteDevice::NaN gives undefined errors when compiled by a regular C++ compiler, hence the earlier guarding requirement 
    inline DEVICE const double NaN = __longlong_as_double(0xfff8000000000000ULL);
#endif
    template<typename T> HOST DEVICE inline constexpr T Const(double v) { return static_cast<T>(v); }

    // 1-Argument Wrapper
    #define CUDA_MATH_1ARG(FuncName, FloatFunc, DoubleFunc) \
    template<typename T> HOST DEVICE FORCE_INLINE T FuncName(T x) { \
        if constexpr (std::is_same_v<T, float>) return FloatFunc(x); \
        else return DoubleFunc(x); \
    }

    // 2-Argument Wrapper
    #define CUDA_MATH_2ARG(FuncName, FloatFunc, DoubleFunc) \
    template<typename T1, typename T2> \
    HOST DEVICE FORCE_INLINE auto FuncName(T1 a, T2 b) { \
        if constexpr (std::is_same_v<T1, float>) { \
            return FloatFunc(static_cast<float>(a), static_cast<float>(b)); \
        } else { \
            return DoubleFunc(static_cast<double>(a), static_cast<double>(b)); \
        } \
    }

    // Reciprocal Wrapper (1.0 / Function)
    #define CUDA_MATH_RECIP(FuncName, BaseFunc) \
    template<typename T> HOST DEVICE FORCE_INLINE T FuncName(T x) { \
        return T(1.0) / BaseFunc(x); \
    }

    // Inverse Reciprocal Wrapper (Function(1.0 / x))
    #define CUDA_MATH_INV_RECIP(FuncName, BaseFunc) \
    template<typename T> HOST DEVICE FORCE_INLINE T FuncName(T x) { \
        return BaseFunc(T(1.0) / x); \
    }

    CUDA_MATH_1ARG(Abs,   fabsf,  fabs)
    CUDA_MATH_1ARG(Acos,  acosf,  acos)
    CUDA_MATH_1ARG(Acosh, acoshf, acosh)
    CUDA_MATH_1ARG(Asin,  asinf,  asin)
    CUDA_MATH_1ARG(Asinh, asinhf, asinh)
    CUDA_MATH_1ARG(Atan,  atanf,  atan)
    CUDA_MATH_1ARG(Atanh, atanhf, atanh)
    CUDA_MATH_1ARG(Ceil,  ceilf,  ceil)
    CUDA_MATH_1ARG(Cos,   cosf,   cos)
    CUDA_MATH_1ARG(Cosh,  coshf,  cosh)
    CUDA_MATH_1ARG(Exp,   expf,   exp)
    CUDA_MATH_1ARG(Expm1, expm1f, expm1)
    CUDA_MATH_1ARG(Floor, floorf, floor)
    CUDA_MATH_1ARG(Log,   logf,   log)
    CUDA_MATH_1ARG(Log10, log10f, log10)
    CUDA_MATH_1ARG(Log1p, log1pf, log1p)
    CUDA_MATH_1ARG(Log2,  log2f,  log2)
    CUDA_MATH_1ARG(Sin,   sinf,   sin)
    CUDA_MATH_1ARG(Sinh,  sinhf,  sinh)
    CUDA_MATH_1ARG(Sqrt,  sqrtf,  sqrt)
    CUDA_MATH_1ARG(Tan,   tanf,   tan)
    CUDA_MATH_1ARG(Tanh,  tanhf,  tanh)

    CUDA_MATH_2ARG(Atan2, atan2f, atan2)
    CUDA_MATH_2ARG(Pow,   powf,   pow)

    CUDA_MATH_RECIP(Sec,  Cos)
    CUDA_MATH_RECIP(Csc,  Sin)
    CUDA_MATH_RECIP(Cot,  Tan)
    CUDA_MATH_RECIP(Sech, Cosh)
    CUDA_MATH_RECIP(Csch, Sinh)
    CUDA_MATH_RECIP(Coth, Tanh)

    CUDA_MATH_INV_RECIP(Asec,  Acos)
    CUDA_MATH_INV_RECIP(Acsc,  Asin)
    CUDA_MATH_INV_RECIP(Acot,  Atan)
    CUDA_MATH_INV_RECIP(Asech, Acosh)
    CUDA_MATH_INV_RECIP(Acsch, Asinh)
    CUDA_MATH_INV_RECIP(Acoth, Atanh)

    template<typename T> HOST DEVICE FORCE_INLINE T Sign(T x) {
        return static_cast<T>((T(0) < x) - (x < T(0)));
    }

    #undef CUDA_MATH_1ARG
    #undef CUDA_MATH_2ARG
    #undef CUDA_MATH_RECIP
    #undef CUDA_MATH_INV_RECIP
}

#endif /* CHASTEMATHDEVICENAMESPACES_HPP_ */