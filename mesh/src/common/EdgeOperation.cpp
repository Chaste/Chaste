/*

Copyright (c) 2005-2023, University of Oxford.
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

#include "EdgeOperation.hpp"

EdgeOperation::EdgeOperation()
{
}

EdgeOperation::EdgeOperation(EDGE_OPERATION operation,
                             unsigned elementIndex,
                             EdgeRemapInfo remapInfo,
                             const bool isIndexRemapped)
    : mOperation(operation),
      mElementIndex(elementIndex),
      mElementIndex2(0),
      mRemapInfo(std::move(remapInfo)),
      mRemapInfo2(),
      mIsElementIndexRemapped(isIndexRemapped)
{
}

EdgeOperation::EdgeOperation(unsigned elementIndex,
                             unsigned elementIndex2,
                             EdgeRemapInfo remapInfo,
                             EdgeRemapInfo remapInfo2)
    : mOperation(EDGE_OPERATION_DIVIDE),
      mElementIndex(elementIndex),
      mElementIndex2(elementIndex2),
      mRemapInfo(std::move(remapInfo)),
      mRemapInfo2(std::move(remapInfo2)),
      mIsElementIndexRemapped(false)
{
}

EDGE_OPERATION EdgeOperation::GetOperation() const
{
    return mOperation;
}

unsigned EdgeOperation::GetElementIndex() const
{
    return this->mElementIndex;
}

void EdgeOperation::SetElementIndex(const unsigned index)
{
    this->mElementIndex  = index;
}

unsigned EdgeOperation::GetElementIndex2() const
{
    return this->mElementIndex2;
}

void EdgeOperation::SetElementIndex2(const unsigned index)
{
    this->mElementIndex2 = index;
}

const EdgeRemapInfo& EdgeOperation::rGetRemapInfo() const
{
    return mRemapInfo;
}

const EdgeRemapInfo& EdgeOperation::rGetRemapInfo2() const
{
    return mRemapInfo2;
}

bool EdgeOperation::IsElementIndexRemapped() const
{
    return mIsElementIndexRemapped;
}