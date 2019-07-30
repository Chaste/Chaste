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


#ifndef TESTPLANESTIMULUSCELLFACTORY_HPP_
#define TESTPLANESTIMULUSCELLFACTORY_HPP_

#include "TetrahedralMesh.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "HeartGeometryInformation.hpp"
#include "HeartConfig.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestPlaneStimulusCellFactory : public CxxTest::TestSuite
{
public:
    void TestBasicContructor()
    {
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(2,2,2);
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> cell_factory1;
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> cell_factory2(-100); //  stimulus voltage

        //Try getting before setting. It should throw (covers the exception)
        TS_ASSERT_THROWS_THIS(cell_factory1.GetMesh(),"The mesh object has not been set in the cell factory");

        cell_factory1.SetMesh(&mesh);
        cell_factory2.SetMesh(&mesh);

        for (unsigned node_num=0; node_num<mesh.GetNumNodes(); node_num++)
        {
            Node<3>* p_node=mesh.GetNode(node_num);

            HeartConfig::Instance()->SetOdeTimeStep(0.01);
            AbstractCardiacCellInterface* p_cell1=cell_factory1.CreateCardiacCellForTissueNode(p_node);

            HeartConfig::Instance()->SetOdeTimeStep(0.001);
            AbstractCardiacCellInterface* p_cell2=cell_factory2.CreateCardiacCellForTissueNode(p_node);
            // compute 1 second
            p_cell1->Compute(0.0,1.0);
            p_cell2->Compute(0.0,1.0);
            if (p_node->GetPoint()[0] == 0.0)
            {
                // stimulated cell
                TS_ASSERT_DELTA(p_cell1->GetVoltage(), 161.1, 0.1);
                TS_ASSERT_DELTA(p_cell2->GetVoltage(), 42.0, 0.1);
            }
            else
            {
                // unstimulated cell
                TS_ASSERT_DELTA(p_cell1->GetVoltage(), -83.8, 0.1);
                TS_ASSERT_DELTA(p_cell2->GetVoltage(), -83.8, 0.1);
            }
            delete p_cell1;
            delete p_cell2;
        }
    }

    void TestHeartGeometryIntoCellFactory()
    {
        TetrahedralMesh<2,2> mesh;
        //This mesh will have 6 nodes per face, spaced by 1
        mesh.ConstructRectangularMesh(5, 5);

        OutputFileHandler handler("CellFactory", false);
        std::string left_file="left";
        std::string right_file="right";
        if (PetscTools::AmMaster())
        {
            out_stream p_left_file = handler.OpenOutputFile(left_file);
            out_stream p_right_file = handler.OpenOutputFile(right_file);

            for (unsigned index=0; index<mesh.GetNumNodes(); index++)
            {
                // Get the nodes at the left face of the square
                if (fabs(mesh.GetNode(index)->rGetLocation()[0]) < 1e-6)
                {
                    *p_left_file<< index <<"\n";
                }
                // Get the nodes at the right face of the square
                if (fabs(mesh.GetNode(index)->rGetLocation()[0]-5.0) < 1e-6)
                {
                    *p_right_file<< index <<"\n";
                }

            }
        }

        HeartGeometryInformation<2> info(mesh, handler.GetOutputDirectoryFullPath()+left_file, handler.GetOutputDirectoryFullPath()+right_file, true);

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory_1;

        HeartGeometryInformation<2>* p_info_from_cell_factory = NULL;
        //get it from the cell factory before setting it, it should throw...(covers the exception)
        TS_ASSERT_THROWS_THIS(cell_factory_1.GetHeartGeometryInformation(), "HeartGeometryInformation object has not been set in the cell factory");

        //set the heart geometry information
        cell_factory_1.SetHeartGeometryInformation(&info);
        //now we get it from the cell factory
        p_info_from_cell_factory = cell_factory_1.GetHeartGeometryInformation();

        //check that the object obtained from the cell factory is the same as the one created above.
        for (unsigned index=0; index<mesh.GetNumNodes(); index++)
        {
            double x = mesh.GetNode(index)->rGetLocation()[0];
            TS_ASSERT_EQUALS(info.CalculateRelativeWallPosition(index),(5.0-x)/5.0);
            TS_ASSERT_DELTA(info.CalculateRelativeWallPosition(index),p_info_from_cell_factory->CalculateRelativeWallPosition(index), 1e-15);
            TS_ASSERT_EQUALS(info.rGetDistanceMapEpicardium()[index],x);
            TS_ASSERT_EQUALS(info.rGetDistanceMapEpicardium()[index],p_info_from_cell_factory->rGetDistanceMapEpicardium()[index]);
            TS_ASSERT_EQUALS(info.rGetDistanceMapEndocardium()[index],(5.0-x));
            TS_ASSERT_EQUALS(info.rGetDistanceMapEndocardium()[index],p_info_from_cell_factory->rGetDistanceMapEndocardium()[index]);
        }
    }
};

#endif /*TESTPLANESTIMULUSCELLFACTORY_HPP_*/
