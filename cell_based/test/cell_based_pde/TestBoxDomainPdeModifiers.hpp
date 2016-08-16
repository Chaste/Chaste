/*

Copyright (c) 2005-2016, University of Oxford.
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

#ifndef TESTBOXDOMAINPDEMODIFIERS_HPP_
#define TESTBOXDOMAINPDEMODIFIERS_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "CheckpointArchiveTypes.hpp"
#include "EllipticBoxDomainPdeModifier.hpp"
#include "ParabolicBoxDomainPdeModifier.hpp"
#include "PdeAndBoundaryConditions.hpp"
#include "ParabolicPdeAndBoundaryConditions.hpp"
#include "SimpleUniformSourcePde.hpp"
#include "UniformSourceParabolicPde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "ArchiveOpener.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestBoxDomainPdeModifiers : public AbstractCellBasedWithTimingsTestSuite
{
public:
    void TestEllipticConstructor() throw(Exception)
    {
        // Make the PDE and BCs
        SimpleUniformSourcePde<2> pde(-0.1);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_and_bc.SetDependentVariableName("averaged quantity");

        // Make domain
        ChastePoint<2> lower(-10.0, -10.0);
        ChastePoint<2> upper(10.0, 10.0);
        ChasteCuboid<2> cuboid(lower, upper);

        // Create an elliptic PDE modifier object using this PDE and BCs object
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (&pde_and_bc, cuboid, 2.0));

        // Test that member variables are initialised correctly
        TS_ASSERT_EQUALS(p_pde_modifier->mpPdeAndBcs->rGetDependentVariableName(), "averaged quantity");

        // Check mesh
        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumNodes(),121u);
        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumBoundaryNodes(),40u);
        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumElements(),200u);
        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumBoundaryElements(),40u);

        ChasteCuboid<2> bounding_box = p_pde_modifier->mpFeMesh->CalculateBoundingBox();
        TS_ASSERT_DELTA(bounding_box.rGetUpperCorner()[0],10,1e-5);
        TS_ASSERT_DELTA(bounding_box.rGetUpperCorner()[1],10,1e-5);
        TS_ASSERT_DELTA(bounding_box.rGetLowerCorner()[0],-10,1e-5);
        TS_ASSERT_DELTA(bounding_box.rGetLowerCorner()[1],-10,1e-5);

        // Coverage
        TS_ASSERT_EQUALS(p_pde_modifier->GetOutputGradient(),false); // Defaults to false
        p_pde_modifier->SetOutputGradient(true);
        TS_ASSERT_EQUALS(p_pde_modifier->GetOutputGradient(),true);
    }

    void TestParabolicConstructor() throw(Exception)
    {
        // Make the PDE and BCs
        UniformSourceParabolicPde<2> pde(-0.1);
        ConstBoundaryCondition<2> bc(1.0);
        ParabolicPdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_and_bc.SetDependentVariableName("averaged quantity");

        // Make domain
        ChastePoint<2> lower(-10.0, -10.0);
        ChastePoint<2> upper(10.0, 10.0);
        ChasteCuboid<2> cuboid(lower, upper);

        // Create a parabolic PDE modifier object using this PDE and BCs object
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<2>, p_pde_modifier, (&pde_and_bc,cuboid,2.0));

        // Test that member variables are initialised correctly
        TS_ASSERT_EQUALS(p_pde_modifier->mpPdeAndBcs->rGetDependentVariableName(), "averaged quantity");

        // Check mesh
        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumNodes(),121u);
        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumBoundaryNodes(),40u);
        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumElements(),200u);
        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumBoundaryElements(),40u);

        ChasteCuboid<2> bounding_box = p_pde_modifier->mpFeMesh->CalculateBoundingBox();
        TS_ASSERT_DELTA(bounding_box.rGetUpperCorner()[0],10,1e-5);
        TS_ASSERT_DELTA(bounding_box.rGetUpperCorner()[1],10,1e-5);
        TS_ASSERT_DELTA(bounding_box.rGetLowerCorner()[0],-10,1e-5);
        TS_ASSERT_DELTA(bounding_box.rGetLowerCorner()[1],-10,1e-5);

        // Coverage
        TS_ASSERT_EQUALS(p_pde_modifier->GetOutputGradient(),false); // Defaults to false
        p_pde_modifier->SetOutputGradient(true);
        TS_ASSERT_EQUALS(p_pde_modifier->GetOutputGradient(),true);
    }

    void TestArchiveEllipticBoxDomainPdeModifier() throw(Exception)
    {
        // Create a file for archiving
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "EllipticBoxDomainPdeModifier.arch";

        // Separate scope to write the archive
        {
            // Make the PDE and BCs
            SimpleUniformSourcePde<2> pde(-0.1);
            ConstBoundaryCondition<2> bc(1.0);
            PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
            pde_and_bc.SetDependentVariableName("averaged quantity");

            // Make domain
            ChastePoint<2> lower(-10.0, -10.0);
            ChastePoint<2> upper(10.0, 10.0);
            ChasteCuboid<2> cuboid(lower, upper);

            // Initialise an elliptic PDE modifier object using this PDE and BCs object
            AbstractCellBasedSimulationModifier<2,2>* const p_modifier = new EllipticBoxDomainPdeModifier<2>(&pde_and_bc, cuboid, 2.0);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            output_arch << p_modifier;
            delete p_modifier;
        }

        // Separate scope to read the archive
        {
            AbstractCellBasedSimulationModifier<2,2>* p_modifier2;

            // Restore the modifier
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_modifier2;

            // See whether we read out the correct variable name area
            std::string variable_name = (static_cast<EllipticBoxDomainPdeModifier<2>*>(p_modifier2))->mpPdeAndBcs->rGetDependentVariableName();
            TS_ASSERT_EQUALS(variable_name, "averaged quantity");

            ///\todo #2687 The archive created a new PDE object. Memory-management of mpPdeAndBcs is not enabled. Suggest using a shared-pointer.
            delete (static_cast<EllipticBoxDomainPdeModifier<2>*>(p_modifier2))->mpPdeAndBcs;
            delete p_modifier2;
        }
    }

    void TestArchiveParabolicBoxDomainPdeModifier() throw(Exception)
    {
        // Create a file for archiving
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ParabolicBoxDomainPdeModifier.arch";

        // Separate scope to write the archive
        {
            // Make the PDE and BCs
            UniformSourceParabolicPde<2> pde(-0.1);
            ConstBoundaryCondition<2> bc(1.0);
            ParabolicPdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
            pde_and_bc.SetDependentVariableName("averaged quantity");

            ///\todo #2687
            // Make a dummy solution (because the archiver expects to be able to read/write PETSc Vecs).
            Vec vector = PetscTools::CreateAndSetVec(10, -42.0);
            pde_and_bc.SetSolution(vector);

            // Make domain
            ChastePoint<2> lower(-10.0, -10.0);
            ChastePoint<2> upper(10.0, 10.0);
            ChasteCuboid<2> cuboid(lower, upper);

            // Initialise a parabolic PDE modifier object using this PDE and BCs object
            AbstractCellBasedSimulationModifier<2,2>* const p_modifier = new ParabolicBoxDomainPdeModifier<2>(&pde_and_bc,cuboid,2.0);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            output_arch << p_modifier;
            delete p_modifier;
        }

        // Separate scope to read the archive
        {
            AbstractCellBasedSimulationModifier<2,2>* p_modifier2;

            // Restore the modifier
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_modifier2;

            // See whether we read out the correct variable name
            std::string variable_name = (static_cast<ParabolicBoxDomainPdeModifier<2>*>(p_modifier2))->mpPdeAndBcs->rGetDependentVariableName();
            TS_ASSERT_EQUALS(variable_name, "averaged quantity");

            ///\todo #2687 The archive created a new PDE object. Memory-management of mpPdeAndBcs is not enabled. Suggest using a shared-pointer.
            delete (static_cast<ParabolicBoxDomainPdeModifier<2>*>(p_modifier2))->mpPdeAndBcs;
            delete p_modifier2;
        }
    }
};

#endif /*TESTBOXDOMAINPDEMODIFIERS_HPP_*/
