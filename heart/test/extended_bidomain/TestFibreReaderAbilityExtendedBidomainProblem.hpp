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

#include <cxxtest/TestSuite.h>
#include "ExtendedBidomainProblem.hpp"
#include "CorriasBuistSMCModified.hpp"
#include "BidomainProblem.hpp"
#include "PetscSetupAndFinalize.hpp"

#ifndef TESTFIBREREADERABILITYEXTENDEDBIDOMAINPROBLEM_HPP_
#define TESTFIBREREADERABILITYEXTENDEDBIDOMAINPROBLEM_HPP_

class UnStimulated3DCellFactory: public AbstractCardiacCellFactory<3>
{

public:
    UnStimulated3DCellFactory() : AbstractCardiacCellFactory<3>()
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
    {
        CorriasBuistSMCModified *cell;
        cell = new CorriasBuistSMCModified(mpSolver, mpZeroStimulus);

        cell->SetFakeIccStimulusPresent(false);
        return cell;
    }
};

class UnStimulatedCellFactory: public AbstractCardiacCellFactory<2>
{

public:
    UnStimulatedCellFactory() : AbstractCardiacCellFactory<2>()
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        CorriasBuistSMCModified *cell;
        cell = new CorriasBuistSMCModified(mpSolver, mpZeroStimulus);

        cell->SetFakeIccStimulusPresent(false);//it will get it from the real ICC, via gap junction
        return cell;
    }
};

class TestFibreReaderAbilityExtendedBidomainProblem: public CxxTest::TestSuite
{

public:

    void TestFibreAbilityNoFibreExtendedProblem()
    {
        HeartConfig::Instance()->Reset();
        HeartEventHandler::Instance()->Reset();

        UnStimulatedCellFactory un_stimulated_cell_factory;
        ExtendedBidomainProblem<2> extended_problem( &un_stimulated_cell_factory , &un_stimulated_cell_factory);
        extended_problem.SetExtendedBidomainParameters(1, 2, 3, 4, 5, 6);

        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1, 2));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(3, 4));
        extended_problem.SetIntracellularConductivitiesForSecondCell(Create_c_vector(5, 6));

        // This mesh has parallel to x axis fibre for x < 0.5 mm and aligned 45 degrees for x > 0.5 mm. However,
        // in this case the fibre is not read.
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/2D_0_to_1mm_800_elements");

        extended_problem.Initialise();

        ExtendedBidomainTissue<2>* extended_tissue = extended_problem.GetExtendedBidomainTissue();

        c_matrix<double, 2, 2> modified_extracellular;
        c_matrix<double, 2, 2> modified_intracellular;
        c_matrix<double, 2, 2> modified_intracellular_second_cell;
        double tol = 1e-5;

        // We should expect unchanged tensors for all elements.
        for (AbstractTetrahedralMesh<2,2>::ElementIterator it = extended_problem.rGetMesh().GetElementIteratorBegin();
                it != extended_problem.rGetMesh().GetElementIteratorBegin();
                ++it)
        {
            modified_extracellular =  extended_tissue->rGetExtracellularConductivityTensor(it->GetIndex());
            TS_ASSERT_DELTA(modified_extracellular(0,0), 3, tol);
            TS_ASSERT_DELTA(modified_extracellular(0,1), 0, tol);
            TS_ASSERT_DELTA(modified_extracellular(1,0), 0, tol);
            TS_ASSERT_DELTA(modified_extracellular(1,1), 4, tol);

            modified_intracellular =  extended_tissue->rGetIntracellularConductivityTensor(it->GetIndex());
            TS_ASSERT_DELTA(modified_intracellular(0,0), 1, tol);
            TS_ASSERT_DELTA(modified_intracellular(0,1), 0, tol);
            TS_ASSERT_DELTA(modified_intracellular(1,0), 0, tol);
            TS_ASSERT_DELTA(modified_intracellular(1,1), 2, tol);

            modified_intracellular_second_cell =  extended_tissue->rGetIntracellularConductivityTensorSecondCell(it->GetIndex());
            TS_ASSERT_DELTA(modified_intracellular_second_cell(0,0), 5, tol);
            TS_ASSERT_DELTA(modified_intracellular_second_cell(0,1), 0, tol);
            TS_ASSERT_DELTA(modified_intracellular_second_cell(1,0), 0, tol);
            TS_ASSERT_DELTA(modified_intracellular_second_cell(1,1), 6, tol);
        }
    }

    void TestOrthoFibreAbilityExtendedProblem()
    {
        HeartConfig::Instance()->Reset();
        HeartEventHandler::Instance()->Reset();

        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1, 2));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(3, 4));

        // This mesh has parallel to x axis fibre for x < 0.5 mm and aligned 45 degrees for x > 0.5 mm.
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/2D_0_to_1mm_800_elements", cp::media_type::Orthotropic);

        UnStimulatedCellFactory un_stimulated_cell_factory;

        ExtendedBidomainProblem<2> extended_problem( &un_stimulated_cell_factory , &un_stimulated_cell_factory);

        extended_problem.SetExtendedBidomainParameters(1, 2, 3, 4, 5, 6);

        extended_problem.SetIntracellularConductivitiesForSecondCell(Create_c_vector(5, 6));

        extended_problem.Initialise();

        ExtendedBidomainTissue<2>* extended_tissue = extended_problem.GetExtendedBidomainTissue();

        c_matrix<double, 2, 2> modified_extracellular;
        c_matrix<double, 2, 2> modified_intracellular;
        c_matrix<double, 2, 2> modified_intracellular_second_cell;
        double tol = 1e-5;

        // We should expect unchanged tensors for elements with x < 0.5 mm and modified tensors for x > 0.5
        for (AbstractTetrahedralMesh<2,2>::ElementIterator it = extended_problem.rGetMesh().GetElementIteratorBegin();
             it != extended_problem.rGetMesh().GetElementIteratorBegin();
             ++it)
        {
            // Get effective conductivity tensors for element < 0.5 mm
            if (it->CalculateCentroid()[0] < 0.05)
            {
                modified_extracellular =  extended_tissue->rGetExtracellularConductivityTensor(it->GetIndex());
                TS_ASSERT_DELTA(modified_extracellular(0,0), 3, tol);
                TS_ASSERT_DELTA(modified_extracellular(0,1), 0, tol);
                TS_ASSERT_DELTA(modified_extracellular(1,0), 0, tol);
                TS_ASSERT_DELTA(modified_extracellular(1,1), 4, tol);

                modified_intracellular =  extended_tissue->rGetIntracellularConductivityTensor(it->GetIndex());
                TS_ASSERT_DELTA(modified_intracellular(0,0), 1, tol);
                TS_ASSERT_DELTA(modified_intracellular(0,1), 0, tol);
                TS_ASSERT_DELTA(modified_intracellular(1,0), 0, tol);
                TS_ASSERT_DELTA(modified_intracellular(1,1), 2, tol);

                modified_intracellular_second_cell =  extended_tissue->rGetIntracellularConductivityTensorSecondCell(it->GetIndex());
                TS_ASSERT_DELTA(modified_intracellular_second_cell(0,0), 5, tol);
                TS_ASSERT_DELTA(modified_intracellular_second_cell(0,1), 0, tol);
                TS_ASSERT_DELTA(modified_intracellular_second_cell(1,0), 0, tol);
                TS_ASSERT_DELTA(modified_intracellular_second_cell(1,1), 6, tol);
            }
            else
            {
                modified_extracellular =  extended_tissue->rGetExtracellularConductivityTensor(it->GetIndex());
                TS_ASSERT_DELTA(modified_extracellular(0,0), 3.5, tol);
                TS_ASSERT_DELTA(modified_extracellular(0,1), -0.5, tol);
                TS_ASSERT_DELTA(modified_extracellular(1,0), -0.5, tol);
                TS_ASSERT_DELTA(modified_extracellular(1,1), 3.5, tol);

                modified_intracellular =  extended_tissue->rGetIntracellularConductivityTensor(it->GetIndex());
                TS_ASSERT_DELTA(modified_intracellular(0,0), 1.5, tol);
                TS_ASSERT_DELTA(modified_intracellular(0,1), -0.5, tol);
                TS_ASSERT_DELTA(modified_intracellular(1,0), -0.5, tol);
                TS_ASSERT_DELTA(modified_intracellular(1,1), 1.5, tol);

                modified_intracellular_second_cell =  extended_tissue->rGetIntracellularConductivityTensorSecondCell(it->GetIndex());
                TS_ASSERT_DELTA(modified_intracellular_second_cell(0,0), 5.5, tol);
                TS_ASSERT_DELTA(modified_intracellular_second_cell(0,1), -0.5, tol);
                TS_ASSERT_DELTA(modified_intracellular_second_cell(1,0), -0.5, tol);
                TS_ASSERT_DELTA(modified_intracellular_second_cell(1,1), 5.5, tol);
            }
        }
    }

    void Test3DAxiFibreAbilityExtendedProblem()
    {
        HeartConfig::Instance()->Reset();
        HeartEventHandler::Instance()->Reset();
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1, 2, 2));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(3, 4, 4));

        // This mesh defines a single tetrahdera with fibre oriented along y axis
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/3D_Single_tetrahedron_element", cp::media_type::Axisymmetric);

        // We ignore mesh partitioning as this mesh has only one element.
        HeartConfig::Instance()->SetMeshPartitioning("dumb");

        UnStimulated3DCellFactory un_stimulated_3Dcell_factory;

        ExtendedBidomainProblem<3> extended_problem( &un_stimulated_3Dcell_factory , &un_stimulated_3Dcell_factory);
        extended_problem.SetExtendedBidomainParameters(1, 2, 3, 4, 5, 6);
        extended_problem.SetIntracellularConductivitiesForSecondCell(Create_c_vector(5, 6, 6));
        extended_problem.Initialise();

        ExtendedBidomainTissue<3>* extended_tissue = extended_problem.GetExtendedBidomainTissue();

        c_matrix<double, 3, 3> modified_extracellular;
        c_matrix<double, 3, 3> modified_intracellular;
        c_matrix<double, 3, 3> modified_intracellular_second_cell;

        double tol = 1e-5;

        // We should expect modified tensors for the mesh.
        for (AbstractTetrahedralMesh<3,3>::ElementIterator it = extended_problem.rGetMesh().GetElementIteratorBegin();
                it != extended_problem.rGetMesh().GetElementIteratorBegin();
                ++it)
        {
            modified_extracellular =  extended_tissue->rGetExtracellularConductivityTensor(it->GetIndex());
            TS_ASSERT_DELTA(modified_extracellular(0,0), 4, tol);
            TS_ASSERT_DELTA(modified_extracellular(0,1), 0, tol);
            TS_ASSERT_DELTA(modified_extracellular(0,2), 0, tol);
            TS_ASSERT_DELTA(modified_extracellular(1,0), 0, tol);
            TS_ASSERT_DELTA(modified_extracellular(1,1), 3, tol);
            TS_ASSERT_DELTA(modified_extracellular(1,2), 0, tol);
            TS_ASSERT_DELTA(modified_extracellular(2,0), 0, tol);
            TS_ASSERT_DELTA(modified_extracellular(2,1), 0, tol);
            TS_ASSERT_DELTA(modified_extracellular(2,2), 4, tol);

            modified_intracellular =  extended_tissue->rGetIntracellularConductivityTensor(it->GetIndex());
            TS_ASSERT_DELTA(modified_intracellular(0,0), 2, tol);
            TS_ASSERT_DELTA(modified_intracellular(0,1), 0, tol);
            TS_ASSERT_DELTA(modified_intracellular(0,2), 0, tol);
            TS_ASSERT_DELTA(modified_intracellular(1,0), 0, tol);
            TS_ASSERT_DELTA(modified_intracellular(1,1), 1, tol);
            TS_ASSERT_DELTA(modified_intracellular(1,2), 0, tol);
            TS_ASSERT_DELTA(modified_intracellular(2,0), 0, tol);
            TS_ASSERT_DELTA(modified_intracellular(2,1), 0, tol);
            TS_ASSERT_DELTA(modified_intracellular(2,2), 2, tol);

            modified_intracellular_second_cell = extended_tissue->rGetIntracellularConductivityTensorSecondCell(it->GetIndex());
            TS_ASSERT_DELTA(modified_intracellular_second_cell(0,0), 6, tol);
            TS_ASSERT_DELTA(modified_intracellular_second_cell(0,1), 0, tol);
            TS_ASSERT_DELTA(modified_intracellular_second_cell(0,2), 0, tol);
            TS_ASSERT_DELTA(modified_intracellular_second_cell(1,0), 0, tol);
            TS_ASSERT_DELTA(modified_intracellular_second_cell(1,1), 5, tol);
            TS_ASSERT_DELTA(modified_intracellular_second_cell(1,2), 0, tol);
            TS_ASSERT_DELTA(modified_intracellular_second_cell(2,0), 0, tol);
            TS_ASSERT_DELTA(modified_intracellular_second_cell(2,1), 0, tol);
            TS_ASSERT_DELTA(modified_intracellular_second_cell(2,2), 6, tol);
        }
    }
};

#endif /* TESTFIBREREADERABILITYEXTENDEDBIDOMAINPROBLEM_HPP_ */
