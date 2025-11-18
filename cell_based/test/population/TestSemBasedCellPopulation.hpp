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

#ifndef TESTSEMBASEDCELLPOPULATION_HPP_
#define TESTSEMBASEDCELLPOPULATION_HPP_


#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "SemBasedCellPopulation.hpp"
#include "SemMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "NoCellCycleModel.hpp"

#include "AbstractCellBasedTestSuite.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestSemBasedCellPopulation : public AbstractCellBasedTestSuite
{
public:

    void TestConstructorsAndDestructor() const
    {
        ///\todo
    }

    void TestGetMesh() const
    {
        ///\todo both rGetMesh() methods
    }

    void TestGetElement(unsigned elementIndex)
    {
        ///\todo
    }

    void TestGetNumNodes()
    {
        ///\todo
    }

    void TestGetLocationOfCellCentre()
    {
        ///\todo
    }

    void TestGetNode()
    {
        ///\todo
    }

    void TestGetNeighbouringLocationIndices()
    {
        ///\todo
    }

    void TestAddNode()
    {
        ///\todo
    }

    void TestSetNode()
    {
        ///\todo
    }

    void TestGetElementCorrespondingToCell()
    {
        ///\todo
    }

    void TestGetVolumeOfCell()
    {
        ///\todo
    }

    void TestOutputCellPopulationParameters()
    {
        ///\todo
    }
    
    void TestValidate()
    {
        ///\todo
    }
    
    void TestSaveAndLoad()
    {
        ///\todo
    }

    void TestGetAndSetMethods()
    {
        SemMeshGenerator generator;
        generator.GenerateSingleCell({0.0, 0.0}, {0.5, 0.5}, {8, 8});
        auto p_mesh = generator.GetMesh();

        c_vector<double, 4> boxCollectionDomain{};
        boxCollectionDomain[0] = -1.0;
        boxCollectionDomain[1] =  1.0;
        boxCollectionDomain[2] = -1.0;
        boxCollectionDomain[3] =  1.0;

        p_mesh->SetUpBoxCollection(0.1, boxCollectionDomain);

        // Assertions
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 1);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 64);

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        SemBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // mOutputNodeRegionToVtk
        {
            // default value is ture
            TS_ASSERT(cell_population.GetOutputNodeRegionToVtk())
            cell_population.SetOutputNodeRegionToVtk(false);
            TS_ASSERT(!cell_population.GetOutputNodeRegionToVtk());
        }
    }
};

#endif /*TESTSEMBASEDCELLPOPULATION_HPP_*/
