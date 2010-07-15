/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef TESTCRYPTSTATISTICS_HPP_
#define TESTCRYPTSTATISTICS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CryptProjectionStatistics.hpp"
#include "CryptProjectionForce.hpp"
#include "TissueSimulation.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "RadialSloughingCellKiller.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "MeshBasedTissueWithGhostNodes.hpp"
#include "WntConcentration.hpp"
#include "WildTypeCellMutationState.hpp"

class TestCryptProjectionStatistics : public AbstractCellBasedTestSuite
{

public:

    void TestGetSection() throw (Exception)
    {
        // Set up tissue
        TissueConfig* p_params = TissueConfig::Instance();
        p_params->SetWntStemThreshold(0.95);

        double a = 0.2;
        double b = 2.0;
        p_params->SetCryptProjectionParameterA(a);
        p_params->SetCryptProjectionParameterB(b);

        int num_cells_depth = 20;
        int num_cells_width = 20;
        unsigned thickness_of_ghost_layer = 3;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, thickness_of_ghost_layer, false);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        ChasteCuboid<2> bounding_box=p_mesh->CalculateBoundingBox();
        double width_of_mesh = (num_cells_width/(num_cells_width+2.0*thickness_of_ghost_layer))*(bounding_box.GetWidth(0));
        double height_of_mesh = (num_cells_depth/(num_cells_depth+2.0*thickness_of_ghost_layer))*(bounding_box.GetWidth(1));

        p_mesh->Translate(-width_of_mesh/2,-height_of_mesh/2);

        std::vector<TissueCellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned i=0; i<location_indices.size(); i++)
        {
            SimpleWntCellCycleModel* p_model = new SimpleWntCellCycleModel;
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(TRANSIT);

            TissueCellPtr p_cell(new TissueCell(p_state, p_model));

            p_cell->InitialiseCellCycleModel();
            double birth_time = - RandomNumberGenerator::Instance()->ranf()*
                                  (p_params->GetTransitCellG1Duration()
                                    + p_params->GetSG2MDuration());
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Make a tissue
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Set up the Wnt gradient
        WntConcentration<2>::Instance()->SetType(RADIAL);
        WntConcentration<2>::Instance()->SetTissue(crypt);

        CryptProjectionStatistics statistics(crypt);

        std::vector<TissueCellPtr> test_section = statistics.GetCryptSection(M_PI/2.0);

        // Test the cells are correct
        TS_ASSERT_EQUALS(test_section.size(), 10u);

        unsigned expected_indices[10] = {350,377,402,429,454,481,506,533,558,585};

        for (unsigned i=0; i<test_section.size(); i++)
        {
            TS_ASSERT_EQUALS(crypt.GetLocationIndexUsingCell(test_section[i]), expected_indices[i]);
        }

        // Create the force law and pass in to a std::list
        CryptProjectionForce crypt_projection_force;
        std::vector<AbstractForce<2>* > force_collection;
        force_collection.push_back(&crypt_projection_force);

        // Make a tissue simulation
        TissueSimulation<2> crypt_projection_simulator(crypt, force_collection, false, false);

        // Create a radial cell killer and pass it in to the tissue simulation
        c_vector<double,2> centre = zero_vector<double>(2);
        double crypt_radius = pow(TissueConfig::Instance()->GetCryptLength()/a, 1.0/b);

        RadialSloughingCellKiller killer(&crypt, centre, crypt_radius);
        crypt_projection_simulator.AddCellKiller(&killer);

        // Set up the simulation
        crypt_projection_simulator.SetOutputDirectory("CryptProjectionStatistics");
        crypt_projection_simulator.SetEndTime(0.25);
        crypt_projection_simulator.Solve();

        statistics.LabelSPhaseCells();

        std::vector<TissueCellPtr> test_section2 = statistics.GetCryptSection();
        std::vector<bool> labelled_cells = statistics.AreCryptSectionCellsLabelled(test_section2);

        TS_ASSERT_EQUALS(test_section2.size(), labelled_cells.size());

        // Five of these cells are labelled - at nodes 207, 208, 232, 254 and 255
        for (unsigned i=0; i<test_section2.size(); i++)
        {
            unsigned node_index = crypt.GetLocationIndexUsingCell(test_section2[i]);

            if (node_index == 207u || node_index == 208u || node_index == 232u || node_index == 254u || node_index == 255u)
            {
                TS_ASSERT_EQUALS(labelled_cells[i], true);
            }
            else
            {
                TS_ASSERT_EQUALS(labelled_cells[i], false);
            }
        }

        crypt_projection_simulator.SetEndTime(0.3);
        crypt_projection_simulator.Solve();

        // Tidy up
        WntConcentration<2>::Destroy();
    }
};

#endif /*TESTCRYPTSTATISTICS_HPP_*/
