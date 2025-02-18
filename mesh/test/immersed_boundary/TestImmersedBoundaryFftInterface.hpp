/*

Copyright (c) 2005-2025, University of Oxford.
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

#ifndef TESTFLUIDSOURCE_HPP_
#define TESTFLUIDSOURCE_HPP_

// Needed for test framework
#include <cxxtest/TestSuite.h>

// Includes from projects/ImmersedBoundary
#include "ImmersedBoundaryFftInterface.hpp"

#include <complex>

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class TestImmersedBoundaryFftInterface : public CxxTest::TestSuite
{
public:

    void TestConstructor()
    {
        // Create a square test mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
        nodes.push_back(new Node<2>(2, true, 0.1, 0.1));
        nodes.push_back(new Node<2>(3, true, 0.0, 0.1));

        std::vector<ImmersedBoundaryElement<2, 2>*> elems;
        elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));
        elems.push_back(new ImmersedBoundaryElement<2, 2>(1, nodes));
        elems.push_back(new ImmersedBoundaryElement<2, 2>(2, nodes));

        ImmersedBoundaryMesh<2,2>* p_mesh = new ImmersedBoundaryMesh<2, 2>(nodes, elems);

        std::vector<double> in = {1.0, 0.5, 0.0, 0.5};
        std::vector<std::complex<double>> comp(in.size());
        std::vector<double> out(in.size());

        ImmersedBoundaryFftInterface<2> fft(p_mesh, in.data(), comp.data(), out.data(), true);
        TS_ASSERT_EQUALS(fft.mRealDims[0], p_mesh->GetNumGridPtsX());

        delete p_mesh;
    }

    void TestFftExecuteForward()
    {
        {
            // Create a square test mesh
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.1));
            nodes.push_back(new Node<2>(3, true, 0.0, 0.1));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));
            elems.push_back(new ImmersedBoundaryElement<2, 2>(1, nodes));
            elems.push_back(new ImmersedBoundaryElement<2, 2>(2, nodes));

            ImmersedBoundaryMesh<2,2>* p_mesh = new ImmersedBoundaryMesh<2, 2>(nodes, elems);

            std::vector<double> in;
            for (unsigned i = 0; i < 49152; ++i)
            {
                in.push_back(0.0);
            }
            std::vector<std::complex<double>> comp(in.size());
            std::vector<double> out(in.size() * 2);
            ImmersedBoundaryFftInterface<2> fft(p_mesh, in.data(), comp.data(), out.data(), true);

            TS_ASSERT_EQUALS(fft.mRealDims[0], 128);
            TS_ASSERT_EQUALS(fft.mRealDims[1], 128);
            TS_ASSERT_EQUALS(fft.mCompDims[0], 128);
            TS_ASSERT_EQUALS(fft.mCompDims[1], 65);
            TS_ASSERT_EQUALS(fft.mRealSep, 16384);
            TS_ASSERT_EQUALS(fft.mCompSep, 8320);
            TS_ASSERT_EQUALS(fft.mRealStride, 8);
            TS_ASSERT_EQUALS(fft.mCompStride, 16);


            fft.FftExecuteForward();
            for (auto i : out)
            {
                TS_ASSERT_DELTA(i, 0.0, 1e-6);
            }

            delete p_mesh;
        }
        {
            // Create a square test mesh
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.1));
            nodes.push_back(new Node<2>(3, true, 0.0, 0.1));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));
            elems.push_back(new ImmersedBoundaryElement<2, 2>(1, nodes));
            elems.push_back(new ImmersedBoundaryElement<2, 2>(2, nodes));

            ImmersedBoundaryMesh<2,2>* p_mesh = new ImmersedBoundaryMesh<2, 2>(nodes, elems, {}, 4, 4);

            std::vector<double> in;
            for (unsigned i = 0; i < 48; ++i)
            {
                in.push_back(i);
            }
            for (unsigned i = 0; i < 48; i += 3)
            {
                in[i] = 0.5;
            }
            std::vector<std::complex<double>> comp(in.size());
            std::vector<double> out(in.size() * 2);
            ImmersedBoundaryFftInterface<2> fft(p_mesh, in.data(), comp.data(), out.data(), true);

            TS_ASSERT_EQUALS(fft.mRealDims[0], 4);
            TS_ASSERT_EQUALS(fft.mRealDims[1], 4);
            TS_ASSERT_EQUALS(fft.mCompDims[0], 4);
            TS_ASSERT_EQUALS(fft.mCompDims[1], 3);
            TS_ASSERT_EQUALS(fft.mRealSep, 16);
            TS_ASSERT_EQUALS(fft.mCompSep, 12);
            TS_ASSERT_EQUALS(fft.mRealStride, 8);
            TS_ASSERT_EQUALS(fft.mCompStride, 16);


            fft.FftExecuteForward();
            for (auto i : out)
            {
                TS_ASSERT_DELTA(i, 0.0, 1e-6);
            }

            delete p_mesh;
        }
    }

    void TestFftExecuteInverse()
    {
        // Create a square test mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
        nodes.push_back(new Node<2>(2, true, 0.1, 0.1));
        nodes.push_back(new Node<2>(3, true, 0.0, 0.1));

        std::vector<ImmersedBoundaryElement<2, 2>*> elems;
        elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));
        elems.push_back(new ImmersedBoundaryElement<2, 2>(1, nodes));
        elems.push_back(new ImmersedBoundaryElement<2, 2>(2, nodes));

        ImmersedBoundaryMesh<2,2>* p_mesh = new ImmersedBoundaryMesh<2, 2>(nodes, elems);

        std::vector<double> in;
        for (unsigned i = 0; i < 49152; ++i)
        {
            in.push_back(0.0);
        }
        std::vector<std::complex<double>> comp(in.size());
        std::vector<double> out(in.size() * 2);
        ImmersedBoundaryFftInterface<2> fft(p_mesh, in.data(), comp.data(), out.data(), true);

        fft.FftExecuteForward();
        fft.FftExecuteInverse();
        for (auto i : out)
        {
            TS_ASSERT_DELTA(i, 0.0, 1e-6);
        }
        for (auto i : in)
        {
            TS_ASSERT_DELTA(i, 0.0, 1e-6);
        }

        delete p_mesh;
    }
};

#endif /*TESTFLUIDSOURCE_HPP_*/
