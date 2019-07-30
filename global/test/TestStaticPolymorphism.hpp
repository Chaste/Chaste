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

#ifndef TESTSTATICPOLYMORPHISM_HPP_
#define TESTSTATICPOLYMORPHISM_HPP_

#include <cxxtest/TestSuite.h>
#include "PetscSetupAndFinalize.hpp"

#include "Timer.hpp"

// Single class. Method() not virtual.
class SingleClass
{
public:
    int Method(int x)
    {
        return x-1;
    }

    void Run()
    {
        int total = 0;
        for (int i=0; i<1e8; i++)
        {
            total += Method(i);
        }
    }

    virtual ~SingleClass()
    {}
};

class DynamicBaseclass
{
public:
    virtual int Method(int x)
    {
        return x-1;
    }

    void Run()
    {
        int total = 0;
        for (int i=0; i<1e8; i++)
        {
            total += Method(i);
        }
    }

    virtual ~DynamicBaseclass()
    {}
};

class DynamicSubclass : public DynamicBaseclass
{
public:
    virtual int Method(int x)
    {
        return x+1;
    }

    virtual ~DynamicSubclass()
    {}
};

class DynamicSubsubclass : public DynamicSubclass
{
public:
    virtual int Method(int x)
    {
        return x+2;
    }

    virtual ~DynamicSubsubclass()
    {}
};

class DynamicSubsubsubclass : public DynamicSubsubclass
{
public:
    virtual int Method(int x)
    {
        return x+3;
    }

    virtual ~DynamicSubsubsubclass()
    {}
};

template <class Derived>
class StaticBaseclass
{
public:
    int Method(int x)
    {
        return x-1;
    }

    void Run()
    {
        int total = 0;
        for (int i=0; i<1e8; i++)
        {
            total += static_cast<Derived*>(this)->Method(i);
        }
    }

    virtual ~StaticBaseclass()
    {}
};

class StaticSubclass : public StaticBaseclass<StaticSubclass>
{
public:
    int Method(int x)
    {
        return x+1;
    }

    virtual ~StaticSubclass()
    {}
};

class TestStaticPolymorphism : public CxxTest::TestSuite
{
public:
    void TestDynamicPolymorphism()
    {
        SingleClass single_class;
        Timer::Reset();
        single_class.Run();
        std::cout << "Single class: " << Timer::GetElapsedTime() << std::endl;

        DynamicBaseclass baseclass;
        Timer::Reset();
        baseclass.Run();
        std::cout << "Dynamic baseclass: " << Timer::GetElapsedTime() << std::endl;

        DynamicSubclass subclass;
        Timer::Reset();
        subclass.Run();
        std::cout << "Dynamic subclass: " << Timer::GetElapsedTime() << std::endl;

        DynamicSubsubsubclass subsubsubclass;
        Timer::Reset();
        subsubsubclass.Run();
        std::cout << "Dynamic subsubsubclass: " << Timer::GetElapsedTime() << std::endl;

        StaticSubclass static_subclass;
        Timer::Reset();
        static_subclass.Run();
        std::cout << "Static subclass: " << Timer::GetElapsedTime() << std::endl;
    }
};

#endif /*TESTSTATICPOLYMORPHISM_HPP_*/
