/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef DEVICETOOLS_HPP_
#define DEVICETOOLS_HPP_

#include <cstddef>

/**
 * Hides CUDA from code compiled by the normal C++ compiler. Code such as CardiacCellCollection uses
 * these functions to allocate and copy GPU memory without including any CUDA headers; the CUDA calls
 * themselves are in DeviceTools.cu, which only the CUDA compiler sees. CUDA errors are raised as
 * Exceptions.
 */
class DeviceTools
{
public:
    /**
     * Allocate device memory.
     *
     * @param numBytes  the number of bytes to allocate
     * @return a pointer to the device memory, to be released with Free, or NULL if numBytes is zero
     */
    static void* Allocate(size_t numBytes);

    /**
     * Release device memory.
     *
     * @param pDevice  a pointer returned by Allocate
     */
    static void Free(void* pDevice);

    /**
     * Copy from host memory to device memory.
     *
     * @param pDevice  where to copy to on the device
     * @param pHost  where to copy from on the host
     * @param numBytes  the number of bytes to copy
     */
    static void CopyToDevice(void* pDevice, const void* pHost, size_t numBytes);

    /**
     * Copy from device memory to host memory.
     *
     * @param pHost  where to copy to on the host
     * @param pDevice  where to copy from on the device
     * @param numBytes  the number of bytes to copy
     */
    static void CopyToHost(void* pHost, const void* pDevice, size_t numBytes);

    /**
     * Wait for the kernel that has just been launched to finish.
     */
    static void WaitForKernel();
};

#endif // DEVICETOOLS_HPP_
