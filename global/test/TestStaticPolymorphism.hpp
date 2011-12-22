/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef TESTSTATICPOLYMORPHISM_HPP_
#define TESTSTATICPOLYMORPHISM_HPP_

#include <cxxtest/TestSuite.h>
#include "PetscSetupAndFinalize.hpp"

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
        double start_time = MPI_Wtime();
        single_class.Run();
        std::cout << "Single class: " << MPI_Wtime()-start_time << std::endl;

        DynamicBaseclass baseclass;
        start_time = MPI_Wtime();
        baseclass.Run();
        std::cout << "Dynamic baseclass: " << MPI_Wtime()-start_time << std::endl;

        DynamicSubclass subclass;
        start_time = MPI_Wtime();
        subclass.Run();
        std::cout << "Dynamic subclass: " << MPI_Wtime()-start_time << std::endl;

        DynamicSubsubsubclass subsubsubclass;
        start_time = MPI_Wtime();
        subsubsubclass.Run();
        std::cout << "Dynamic subsubsubclass: " << MPI_Wtime()-start_time << std::endl;

        StaticSubclass static_subclass;
        start_time = MPI_Wtime();
        static_subclass.Run();
        std::cout << "Static subclass: " << MPI_Wtime()-start_time << std::endl;
    }
};

#endif /*TESTSTATICPOLYMORPHISM_HPP_*/
