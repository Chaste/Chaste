/*

Copyright (C) University of Oxford, 2005-2012

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


#ifndef TESTPLANESTIMULUSCELLFACTORY_HPP_
#define TESTPLANESTIMULUSCELLFACTORY_HPP_

#include "TetrahedralMesh.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "HeartGeometryInformation.hpp"
#include "HeartConfig.hpp"

class TestPlaneStimulusCellFactory : public CxxTest::TestSuite
{
public:
    void TestBasicContructor() throw (Exception)
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
            AbstractCardiacCell* p_cell1=cell_factory1.CreateCardiacCellForTissueNode(node_num);

            HeartConfig::Instance()->SetOdeTimeStep(0.001);
            AbstractCardiacCell* p_cell2=cell_factory2.CreateCardiacCellForTissueNode(node_num);
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
    void TestHeartGeometryIntoCellFactory() throw(Exception)
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
            TS_ASSERT_EQUALS(info.CalculateRelativeWallPosition(index),p_info_from_cell_factory->CalculateRelativeWallPosition(index));
            TS_ASSERT_EQUALS(info.rGetDistanceMapEpicardium()[index],x);
            TS_ASSERT_EQUALS(info.rGetDistanceMapEpicardium()[index],p_info_from_cell_factory->rGetDistanceMapEpicardium()[index]);
            TS_ASSERT_EQUALS(info.rGetDistanceMapEndocardium()[index],(5.0-x));
            TS_ASSERT_EQUALS(info.rGetDistanceMapEndocardium()[index],p_info_from_cell_factory->rGetDistanceMapEndocardium()[index]);
        }
    }
};

#endif /*TESTPLANESTIMULUSCELLFACTORY_HPP_*/
