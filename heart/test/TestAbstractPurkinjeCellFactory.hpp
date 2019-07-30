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

#ifndef TESTABSTRACTPURKINJECELLFACTORY_HPP_
#define TESTABSTRACTPURKINJECELLFACTORY_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractPurkinjeCellFactory.hpp"
#include "TrianglesMeshReader.hpp"
#include "MixedDimensionMesh.hpp"
#include "TetrahedralMesh.hpp"
#include "FakeBathCell.hpp"
#include "SimpleStimulus.hpp"
#include "LuoRudy1991.hpp"
#include "DiFrancescoNoble1985.hpp"
#include "PetscSetupAndFinalize.hpp"

class PurkinjeCellFactory : public AbstractPurkinjeCellFactory<2>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    PurkinjeCellFactory()
        : AbstractPurkinjeCellFactory<2>(),
          mpStimulus(new SimpleStimulus(-6000.0, 0.5))
    {
    }

    AbstractCardiacCellInterface* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        ChastePoint<2> location = pNode->GetPoint();

        if (fabs(location[0])<1e-6)
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }

    AbstractCardiacCellInterface* CreatePurkinjeCellForTissueNode(Node<2>* pNode,
                                                                  AbstractCardiacCellInterface* pCardiacCell)
    {
        return new CellDiFrancescoNoble1985FromCellML(mpSolver, mpZeroStimulus);
    }
};


class TestAbstractPurkinjeCellFactory : public CxxTest::TestSuite
{
public:
    void TestPurkinjeCellFactory()
    {
        TrianglesMeshReader<2,2> reader("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements");
        MixedDimensionMesh<2,2> mixed_mesh;
        mixed_mesh.ConstructFromMeshReader(reader);

        PurkinjeCellFactory cell_factory;
        TS_ASSERT_THROWS_THIS(cell_factory.GetMixedDimensionMesh(),
                              "The mixed dimension mesh object has not been set in the cell factory");
        cell_factory.SetMesh(&mixed_mesh);

        for (AbstractTetrahedralMesh<2,2>::NodeIterator current_node = mixed_mesh.GetNodeIteratorBegin();
             current_node != mixed_mesh.GetNodeIteratorEnd();
             ++current_node)
        {
            AbstractCardiacCellInterface* p_cell = cell_factory.CreatePurkinjeCellForNode( &(*current_node) , NULL);
            double y = current_node->rGetLocation()[1];

            // cable nodes are on y=0.05 (we don't test by index because indices may be permuted in parallel).
            if (fabs(y-0.05) < 1e-8)
            {
                TS_ASSERT(dynamic_cast<CellDiFrancescoNoble1985FromCellML*>(p_cell) != NULL);
            }
            else
            {
                TS_ASSERT(dynamic_cast<FakeBathCell*>(p_cell) != NULL);
            }
            delete p_cell;
        }

        TS_ASSERT_EQUALS(cell_factory.GetMixedDimensionMesh(), &mixed_mesh);
        TrianglesMeshReader<2,2> reader2("mesh/test/data/2D_0_to_1mm_200_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader2);
        TS_ASSERT_THROWS_THIS(cell_factory.SetMesh(&mesh),
                              "AbstractPurkinjeCellFactory must take a MixedDimensionMesh");

    }
};

#endif // TESTABSTRACTPURKINJECELLFACTORY_HPP_
