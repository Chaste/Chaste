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


#ifndef _TESTNODEATTRIBUTES_HPP_
#define _TESTNODEATTRIBUTES_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "OutputFileHandler.hpp"
#include "NodeAttributes.hpp"
#include "PetscTools.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestNodeAttributes : public CxxTest::TestSuite
{
private:
    template<unsigned DIM>
    void TestConstructNodeAttributes(double tol)
    {
        NodeAttributes<DIM> node_attributes;

        TS_ASSERT_EQUALS(node_attributes.rGetAttributes().size(), 0u);
        TS_ASSERT_EQUALS(node_attributes.GetRegion(), 0u);
        TS_ASSERT_DELTA(node_attributes.GetRadius(), 0.0, 1e-4);
        TS_ASSERT_EQUALS(node_attributes.IsParticle(), false);

        for (unsigned i = 0; i < DIM; i++)
        {
            TS_ASSERT_DELTA(node_attributes.rGetAppliedForce()[i], 0.0, tol);
        }
    }

public:

    void TestConstructor1d2d3d()
    {
        TestConstructNodeAttributes<1>(1e-4);
        TestConstructNodeAttributes<2>(1e-4);
        TestConstructNodeAttributes<3>(1e-4);
    }

    void TestAttributesContainer()
    {
        NodeAttributes<3> node_attributes;

        node_attributes.AddAttribute(1.0);
        TS_ASSERT_EQUALS(node_attributes.rGetAttributes().size(), 1u);
        TS_ASSERT_DELTA(node_attributes.rGetAttributes()[0], 1.0, 1e-4);

        node_attributes.AddAttribute(2.0);
        TS_ASSERT_EQUALS(node_attributes.rGetAttributes().size(), 2u);
        TS_ASSERT_DELTA(node_attributes.rGetAttributes()[0], 1.0, 1e-4);
        TS_ASSERT_DELTA(node_attributes.rGetAttributes()[1], 2.0, 1e-4);
    }

    void TestRegion()
    {
        NodeAttributes<3> node_attributes;

        TS_ASSERT_EQUALS(node_attributes.GetRegion(), 0u);

        node_attributes.SetRegion(1);
        TS_ASSERT_EQUALS(node_attributes.GetRegion(), 1u);

        node_attributes.SetRegion(2);
        TS_ASSERT_EQUALS(node_attributes.GetRegion(), 2u);
    }

    void TestAppliedForce()
    {
        NodeAttributes<3> node_attributes;

        TS_ASSERT_DELTA(node_attributes.rGetAppliedForce()[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(node_attributes.rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(node_attributes.rGetAppliedForce()[2], 0.0, 1e-4);

        c_vector<double, 3> force_contribution;
        force_contribution[0] = 1.0;
        force_contribution[1] = -23.6;
        force_contribution[2] = 12345.0;

        node_attributes.AddAppliedForceContribution(force_contribution);

        TS_ASSERT_DELTA(node_attributes.rGetAppliedForce()[0], 1.0, 1e-4);
        TS_ASSERT_DELTA(node_attributes.rGetAppliedForce()[1], -23.6, 1e-4);
        TS_ASSERT_DELTA(node_attributes.rGetAppliedForce()[2], 12345.0, 1e-4);

        node_attributes.AddAppliedForceContribution(force_contribution);

        TS_ASSERT_DELTA(node_attributes.rGetAppliedForce()[0], 2.0, 1e-4);
        TS_ASSERT_DELTA(node_attributes.rGetAppliedForce()[1], -47.2, 1e-4);
        TS_ASSERT_DELTA(node_attributes.rGetAppliedForce()[2], 24690.0, 1e-4);

        node_attributes.ClearAppliedForce();

        TS_ASSERT_DELTA(node_attributes.rGetAppliedForce()[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(node_attributes.rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(node_attributes.rGetAppliedForce()[2], 0.0, 1e-4);
    }

    void TestIsParticle()
    {
        NodeAttributes<2> node_attributes;

        node_attributes.SetIsParticle(true);

        TS_ASSERT(node_attributes.IsParticle());

        node_attributes.SetIsParticle(false);

        TS_ASSERT(!node_attributes.IsParticle());
    }

    void TestRadius()
    {
        NodeAttributes<3> node_attributes;

        TS_ASSERT_DELTA(node_attributes.GetRadius(), 0.0, 1e-4);

        TS_ASSERT_THROWS_THIS(node_attributes.SetRadius(-1.0), "Trying to set node attributes mRadius to a negative value.");

        node_attributes.SetRadius(2.6);
        TS_ASSERT_DELTA(node_attributes.GetRadius(), 2.6, 1e-4);
    }

    void TestArchiving()
    {
        EXIT_IF_PARALLEL;

        OutputFileHandler handler("TestNodeAttributes", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "node_attributes.arch";

        {
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            NodeAttributes<3>* p_node_attributes = new NodeAttributes<3>();
            p_node_attributes->SetRegion(2);
            p_node_attributes->SetRadius(1.24);

            c_vector<double, 3> force_contribution = scalar_vector<double>(3, 1.0);
            p_node_attributes->AddAppliedForceContribution(force_contribution);

            p_node_attributes->SetIsParticle(true);

            p_node_attributes->AddAttribute(6.3);
            p_node_attributes->AddAttribute(4.1);

            NodeAttributes<3>* const p_const_node_attributes = p_node_attributes;

            output_arch << p_const_node_attributes;

            delete p_node_attributes;
        }

        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            NodeAttributes<3>* p_node_attributes;
            input_arch >> p_node_attributes;

            TS_ASSERT_EQUALS(p_node_attributes->GetRegion(), 2u);
            TS_ASSERT_DELTA(p_node_attributes->GetRadius(), 1.24, 1e-4);

            TS_ASSERT_EQUALS(p_node_attributes->rGetAppliedForce().size(), 3u);
            TS_ASSERT_DELTA(p_node_attributes->rGetAppliedForce()[0], 1.0, 1e-4);
            TS_ASSERT_DELTA(p_node_attributes->rGetAppliedForce()[1], 1.0, 1e-4);
            TS_ASSERT_DELTA(p_node_attributes->rGetAppliedForce()[2], 1.0, 1e-4);

            TS_ASSERT(p_node_attributes->IsParticle());

            TS_ASSERT_EQUALS(p_node_attributes->rGetAttributes().size(), 2u);
            TS_ASSERT_DELTA(p_node_attributes->rGetAttributes()[0], 6.3, 1e-4);
            TS_ASSERT_DELTA(p_node_attributes->rGetAttributes()[1], 4.1, 1e-4);

            delete p_node_attributes;
        }
    }

    void TestArchivingNullPointer()
    {
        EXIT_IF_PARALLEL;

        OutputFileHandler handler("TestNodeAttributes", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "null_node_attributes.arch";

        {
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            NodeAttributes<3>* const p_node_attributes = NULL;

            output_arch << p_node_attributes;

            delete p_node_attributes;
        }

        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            NodeAttributes<3>* p_node_attributes;
            input_arch >> p_node_attributes;

            TS_ASSERT(!p_node_attributes);

            delete p_node_attributes;
        }
    }
};

#endif //_TESTNODEATTRIBUTES_HPP_
