/*

Copyright (c) 2005-2019, University of Oxford.
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

#include <cassert>
#include <iostream>
#include "ForTestArchiving.hpp"

BaseClass::BaseClass()
    : mTagInBaseClass(4)
{
}
BaseClass::~BaseClass()
{
}

ChildClass::ChildClass() : mTag(1), mpParent(NULL)
{
    assert(mTagInBaseClass == 4);
}
void ChildClass::SetParent(ParentClass* pParent)
{
    mpParent = pParent;
}

void ChildClass::Hello()
{
    std::cout << "Hello!\n";
}

SubChildClass::SubChildClass() : mSubTag(2)
{
    assert(mTag == 1);
}

void SubChildClass::Hello()
{
    std::cout << "Hello from subclass tagged " << mSubTag << std::endl;
}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(ChildClass)
CHASTE_CLASS_EXPORT(SubChildClass)

ParentClass::ParentClass(ChildClass* pChild)
    : mTag(0), mpChild(pChild)
{
    mpChild->SetParent(this);
}
