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

#ifndef HOSTDEVICEMACROS_HPP_
#define HOSTDEVICEMACROS_HPP_

/*
	Compilation of Chaste happens in multiple stages using multiple compilers
    (when the device features have actaully been merged into the codebase).

	The C++ compiler (eg g++ on Ubuntu) handles the majority of compilation as usual
	in a C++ project. 

    When nvcc (the NVIDIA compiler) is active, __CUDACC__ is defined, this allows
    us to hide syntax that only nvcc understands with a macro
    Similarly when hipcc (the AMD compiler) is active, __HIPCC__ is defined, and
    conveniently, hipcc understands most of the same syntax as nvcc.

    Now both nvcc and hipcc are not actually compilers, but compiler orchestration
    tools, in that they use multiple compilers to compile device code. They compile
    in multiple phases:
    
    1) A host code compilation phase, where they offload the work to a standard C++
    compiler

    2) A device code compilation phase, where they use their own propriatary compilers
    to compile code for GPU. 

    When nvcc is in its "device compilation" phase, __CUDA_ARCH__ is defined. 
    The hipcc "equivalent" is: __HIP_DEVICE_COMPILE__

    Note that there is a slight difference between these, __HIP_DEVICE_COMPILE__ gives
    no information on the architecture we are compiling for.

    These variables allow us to write source code for nvcc/hipcc that changes between
    passes (these are currently unused, and untested in Chaste but this is a common use
    for these variables)

    https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/#gpu-compilation
    https://rocm.docs.amd.com/projects/HIP/en/latest/understand/compilers.html
*/

#if defined(__CUDACC__) || defined(__HIPCC__)
    #define USING_DEVICE_COMPILER 1
#else
    #define USING_DEVICE_COMPILER 0
#endif

#if USING_DEVICE_COMPILER
	#define HOST __host__
	#define DEVICE __device__
	#define GLOBAL __global__
	#define FORCE_INLINE __forceinline__

    // Only defined during device code compilation phase of nvcc/hipcc run
	#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
		#define IN_DEVICE_PASS 1
	#else
		#define IN_DEVICE_PASS 0
	#endif

#else
	#define HOST
	#define DEVICE
	#define GLOBAL
	#if defined(__GNUC__) || defined(__clang__)
        #define FORCE_INLINE inline __attribute__((always_inline))
    #elif defined(_MSC_VER)
        #define FORCE_INLINE __forceinline
    #else
        #define FORCE_INLINE inline
    #endif
#endif

#endif // HOSTDEVICEMACROS_HPP_