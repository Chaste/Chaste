/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef TESTELECTROMECHANICSPROBLEMDEFINITION_HPP_
#define TESTELECTROMECHANICSPROBLEMDEFINITION_HPP_

#include <cxxtest/TestSuite.h>
#include "ElectroMechanicsProblemDefinition.hpp"

class TestElectroMechanicsProblemDefinition : public CxxTest::TestSuite
{
public:
    void TestInterface() throw (Exception)
    {
        QuadraticMesh<2> mesh(0.1, 1.0, 1.0);

        ElectroMechanicsProblemDefinition<2> problem_defn(mesh);

        std::vector<unsigned> fixed_nodes;
        fixed_nodes.push_back(0);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetUseDefaultCardiacMaterialLaw(INCOMPRESSIBLE);


        TS_ASSERT_THROWS_THIS(problem_defn.Validate(), "Timestep for mechanics solve hasn't been set yet");

        problem_defn.SetMechanicsSolveTimestep(1.0);
        TS_ASSERT_EQUALS(problem_defn.GetMechanicsSolveTimestep(), 1.0);

        TS_ASSERT_THROWS_THIS(problem_defn.Validate(), "Contraction model hasn't been set yet");

        problem_defn.SetContractionModel(NASH2004, 0.01);
        TS_ASSERT_EQUALS(problem_defn.GetContractionModel(), NASH2004);
        TS_ASSERT_EQUALS(problem_defn.GetContractionModelOdeTimestep(), 0.01);

        problem_defn.SetUseDefaultCardiacMaterialLaw(INCOMPRESSIBLE);
        NashHunterPoleZeroLaw<2>* p_law = dynamic_cast<NashHunterPoleZeroLaw<2>*>(problem_defn.GetIncompressibleMaterialLaw(0));
        TS_ASSERT(p_law != NULL);

        problem_defn.SetUseDefaultCardiacMaterialLaw(COMPRESSIBLE);
        CompressibleExponentialLaw<2>* p_law_2 = dynamic_cast<CompressibleExponentialLaw<2>*>(problem_defn.GetCompressibleMaterialLaw(0));
        TS_ASSERT(p_law_2 != NULL);

        problem_defn.SetDeformationAffectsElectrophysiology(true,false);
        TS_ASSERT_EQUALS(problem_defn.GetDeformationAffectsConductivity(), true);
        TS_ASSERT_EQUALS(problem_defn.GetDeformationAffectsCellModels(), false);

        TS_ASSERT_THROWS_THIS(problem_defn.Validate(),"Deformation affecting the conductivity is currently not implemented fully for compressible problems");

        problem_defn.SetDeformationAffectsElectrophysiology(false,true);
        TS_ASSERT_EQUALS(problem_defn.GetDeformationAffectsConductivity(), false);
        TS_ASSERT_EQUALS(problem_defn.GetDeformationAffectsCellModels(), true);

        TS_ASSERT_EQUALS(problem_defn.ReadFibreSheetDirectionsFromFile(), false);
        FileFinder ortho_finder("some_file.ortho",RelativeTo::ChasteSourceRoot);
        problem_defn.SetVariableFibreSheetDirectionsFile(ortho_finder, false);
        TS_ASSERT_EQUALS(problem_defn.ReadFibreSheetDirectionsFromFile(), true);
        TS_ASSERT_EQUALS(problem_defn.GetFibreSheetDirectionsFile().GetLeafName(), "some_file.ortho");
        TS_ASSERT_EQUALS(problem_defn.GetFibreSheetDirectionsDefinedPerQuadraturePoint(), false);

        FileFinder orthoquad_finder("some_file.orthoquad",RelativeTo::ChasteSourceRoot);
        problem_defn.SetVariableFibreSheetDirectionsFile(orthoquad_finder, true);
        TS_ASSERT_EQUALS(problem_defn.ReadFibreSheetDirectionsFromFile(), true);
        TS_ASSERT_EQUALS(problem_defn.GetFibreSheetDirectionsFile().GetLeafName(), "some_file.orthoquad");
        TS_ASSERT_EQUALS(problem_defn.GetFibreSheetDirectionsDefinedPerQuadraturePoint(), true);

        TS_ASSERT_EQUALS(problem_defn.GetNumIncrementsForInitialDeformation(), 1u);
        problem_defn.SetNumIncrementsForInitialDeformation(4);
        TS_ASSERT_EQUALS(problem_defn.GetNumIncrementsForInitialDeformation(), 4u);

        TS_ASSERT_THROWS_THIS(problem_defn.SetNumIncrementsForInitialDeformation(0), "Number of increments for initial deformation must be 1 or more");

        // shouldn't throw
        problem_defn.SetDeformationAffectsElectrophysiology(false,false);
        problem_defn.Validate();

        problem_defn.SetDeformationAffectsElectrophysiology(false,true);
        TS_ASSERT_THROWS_CONTAINS(problem_defn.Validate(), "Deformation affecting cell models cannot be done when fibres-sheet");
    }
};

#endif // TESTELECTROMECHANICSPROBLEMDEFINITION_HPP_
