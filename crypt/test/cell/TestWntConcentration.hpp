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

#ifndef TESTWNTCONCENTRATION_HPP_
#define TESTWNTCONCENTRATION_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "ArchiveOpener.hpp"
#include "CellsGenerator.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "WntCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "TrianglesMeshReader.hpp"
#include "WildTypeCellMutationState.hpp"
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/**
 * Note that all these tests call setUp() and tearDown() before running,
 * so if you copy them into a new test suite be sure to copy these methods
 * too.
 */
class TestWntConcentration : public AbstractCellBasedTestSuite
{
public:

    void TestNoWnt()
    {
        WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), false);
        p_wnt->SetType(NONE);
        p_wnt->SetCryptLength(22.0);
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), false);   // NONE does not register as a set up Wnt Gradient (so stem cells are not moved)

        TS_ASSERT_EQUALS(p_wnt->GetType(), NONE);

        double height = 5;
        double wnt_level = 0.0;
        wnt_level = p_wnt->GetWntLevel(height);

        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);

        c_vector<double,2> location;
        location[0] = 1.5;
        location[1] = 2.3;

        TS_ASSERT_DELTA(p_wnt->GetWntGradient(location)[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(p_wnt->GetWntGradient(location)[1], 0.0, 1e-12);

        WntConcentration<2>::Destroy();
    }

    void TestLinearWntConcentration()
    {
        WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), false);
        p_wnt->SetType(LINEAR);
        double crypt_length = 22.0;
        p_wnt->SetCryptLength(crypt_length);

        double height = 100;
        double wnt_level = 0.0;
        wnt_level = p_wnt->GetWntLevel(height);

        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);

        height = -1e-12;    // for cells very close to 0 on negative side.
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 1.0, 1e-9);

        height = 21.0;
        wnt_level = p_wnt->GetWntLevel(height);

        TS_ASSERT_DELTA(wnt_level, 1.0-height/crypt_length, 1e-9);

        // Test GetWntGradient() method
        c_vector<double,2> location;
        location[0] = 1.5;
        location[1] = 2.3;

        TS_ASSERT_DELTA(p_wnt->GetWntGradient(location)[0], 0.0, 1e-12);
        // This should be equal to -1/22 = -0.0454
        TS_ASSERT_DELTA(p_wnt->GetWntGradient(location)[1], -0.0454, 1e-4);

        WntConcentration<2>::Destroy();
    }

    void TestExponentialWntConcentration()
    {
        WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), false);
        p_wnt->SetType(EXPONENTIAL);
        double crypt_length = 22.0;
        p_wnt->SetCryptLength(crypt_length);

        TS_ASSERT_DELTA(p_wnt->GetWntConcentrationParameter(),1.0, 1e-9);

        double height = 100;
        double wnt_level = 0.0;
        wnt_level = p_wnt->GetWntLevel(height);

        // For heights above the top of the crypt (no Wnt)
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);

        height = -1e-12;    // for cells very close to 0 on negative side (very strong Wnt = 1.0);
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 1.0, 1e-9);

        // For normal 'in range' Wnt height
        height = 21.0;
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, exp(-height/crypt_length), 1e-9);

        // For a change in lambda
        p_wnt->SetWntConcentrationParameter(0.5);
        TS_ASSERT_DELTA(p_wnt->GetWntConcentrationParameter(),0.5, 1e-9);
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, exp(-(height/crypt_length)/0.5), 1e-9);

        // Test GetWntGradient() method

        c_vector<double,2> location;
        location[0] = 1.5;
        location[1] = 2.3;

        TS_ASSERT_THROWS_THIS(p_wnt->GetWntGradient(location)[0],
                              "No method to calculate gradient of this Wnt type");

        WntConcentration<2>::Destroy();
    }

    void TestOffsetLinearWntConcentration()
    {
        WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();
        p_wnt->SetType(LINEAR);
        double crypt_length = 22.0;
        p_wnt->SetCryptLength(crypt_length);
        p_wnt->SetWntConcentrationParameter(1.0/3.0);
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), false);

        double height = 100;
        double wnt_level = 0.0;
        wnt_level = p_wnt->GetWntLevel(height);

        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);

        height = -1e-12;
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 1.0, 1e-9);

        height = 21.0;
        wnt_level = p_wnt->GetWntLevel(height);

        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);

        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);

        // under a third of the way up the crypt.
        height = 7.0;
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 1.0 - height/((1.0/3.0)*crypt_length), 1e-9);

        // more than a third of the way up the crypt.
        height = 10.0;
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);

        WntConcentration<2>::Destroy();
    }

    void TestRadialWntConcentration()
    {
        WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();
        p_wnt->SetType(RADIAL);
        double crypt_length = 22.0;
        p_wnt->SetCryptLength(crypt_length);
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), false);   // only fully set up when a cell population is assigned

        // Test GetWntLevel(double) method
        double height = 100;
        double wnt_level = 0.0;

        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-4);

        height = -1e-12;
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 1.0, 1e-4);

        height = 21.0;
        wnt_level = p_wnt->GetWntLevel(height);

        TS_ASSERT_DELTA(wnt_level, 0.0454, 1e-4);

        height = 7.0;
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 0.6818, 1e-4);

        // Test GetWntLevel(CellPtr) method

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Translate mesh so that its centre is at (0,0)
        mesh.Translate(-0.5,-0.5);

        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        boost::shared_ptr<AbstractCellProperty> p_stem_type(new StemCellProliferativeType);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            WntCellCycleModel* p_model = new WntCellCycleModel();
            p_model->SetDimension(2);
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            double birth_time = 0.0 - i;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create a crypt
        MeshBasedCellPopulation<2> crypt(mesh, cells);
        p_wnt->SetCellPopulation(crypt);
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), true);    // fully set up now

        WntConcentration<2>::Destroy();

        WntConcentration<2>::Instance()->SetType(NONE);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        crypt_length = 1.0;
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);
        TS_ASSERT_EQUALS(WntConcentration<2>::Instance()->IsWntSetUp(), false);    // not fully set up now it is a NONE type

        WntConcentration<2>::Destroy();
        WntConcentration<2>::Instance()->SetType(RADIAL);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);

        p_wnt = WntConcentration<2>::Instance();
        p_wnt->SetCryptLength(crypt_length);
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), true);    // set up again

        AbstractCellPopulation<2>::Iterator cell_iter = crypt.Begin();

        double wnt_at_cell0 = p_wnt->GetWntLevel(*cell_iter);

        double a = p_wnt->GetCryptProjectionParameterA();
        double b = p_wnt->GetCryptProjectionParameterB();
        TS_ASSERT_DELTA(a, 0.5, 1e-12);
        TS_ASSERT_DELTA(b, 2.0, 1e-12);

        while (cell_iter != crypt.End())
        {
            TS_ASSERT_DELTA(p_wnt->GetWntLevel(*cell_iter), wnt_at_cell0, 1e-12);

            // Test GetWntGradient(CellPtr) method
            c_vector<double,2> cell_location = crypt.GetLocationOfCellCentre(*cell_iter);
            double r = norm_2(cell_location);

            c_vector<double,2> expected_wnt_gradient;
            expected_wnt_gradient[0] = -cell_location[0]*pow(r,b-1.0)/(a*r);
            expected_wnt_gradient[1] = -cell_location[1]*pow(r,b-1.0)/(a*r);

            TS_ASSERT_DELTA(p_wnt->GetWntGradient(*cell_iter)[0],expected_wnt_gradient[0],1e-6);
            TS_ASSERT_DELTA(p_wnt->GetWntGradient(*cell_iter)[1],expected_wnt_gradient[1],1e-6);

            ++cell_iter;
        }

        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), true);

        WntConcentration<2>::Destroy();
    }

    void TestArchiveWntConcentration()
    {
        // Set up simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Give each a birth time of -node_index, so the age = node_index
        std::vector<CellPtr> cells;
        CellsGenerator<WntCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Work out where to put the archive
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "wnt_concentration.arch";
        ArchiveLocationInfo::SetMeshFilename("wnt_concentration_mesh");

        // Create an output archive
        {
            WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();
            p_wnt->SetType(LINEAR);
            p_wnt->SetCryptLength(22.0);
            p_wnt->SetCryptProjectionParameterA(3.3);
            p_wnt->SetCryptProjectionParameterB(4.4);
            p_wnt->SetCellPopulation(cell_population);

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            SerializableSingleton<WntConcentration<2> >* const p_wrapper = p_wnt->GetSerializationWrapper();
            (*p_arch) << p_wrapper;

            WntConcentration<2>::Destroy();
        }

        {
            WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // Write to the archive
            SerializableSingleton<WntConcentration<2> >* p_wrapper;
            (*p_arch) >> p_wrapper;

            TS_ASSERT_EQUALS(WntConcentration<2>::Instance(), p_wnt);
            TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), true);
            TS_ASSERT_DELTA(p_wnt->GetWntLevel(21.0), 1.0 - 21.0/p_wnt->GetCryptLength(), 1e-9);
            TS_ASSERT_DELTA(p_wnt->GetCryptLength(), 22.0, 1e-12);
            TS_ASSERT_DELTA(p_wnt->GetCryptProjectionParameterA(), 3.3, 1e-12);
            TS_ASSERT_DELTA(p_wnt->GetCryptProjectionParameterB(), 4.4, 1e-12);

            AbstractCellPopulation<2>& arch_cell_population = p_wnt->rGetCellPopulation();
            delete (&arch_cell_population);
        }

        WntConcentration<2>::Destroy();
    }

    void TestSingletonnessOfWntConcentration()
    {
        WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();
        p_wnt->SetType(NONE);
        double crypt_length = 22.0;
        p_wnt->SetCryptLength(crypt_length);

        double height = 5;
        double wnt_level = 0.0;
        wnt_level = p_wnt->GetWntLevel(height);

        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);

        TS_ASSERT_THROWS_THIS(p_wnt->SetType(NONE),"Destroy has not been called");
        TS_ASSERT_THROWS_THIS(p_wnt->SetCryptLength(10.0),"Destroy has not been called");
        WntConcentration<2>::Destroy();

        p_wnt = WntConcentration<2>::Instance();
        p_wnt->SetType(LINEAR);
        p_wnt->SetCryptLength(crypt_length);

        height = 100;
        wnt_level = 0.0;
        wnt_level = p_wnt->GetWntLevel(height);

        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);

        height = -1e-12;    // for cells very close to 0 on negative side.
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 1.0, 1e-9);

        height = 21.0;
        wnt_level = p_wnt->GetWntLevel(height);

        TS_ASSERT_DELTA(wnt_level, 1.0-height/crypt_length, 1e-9);

        TS_ASSERT_THROWS_THIS(p_wnt->SetConstantWntValueForTesting(-10),"WntConcentration<DIM>::SetConstantWntValueForTesting - Wnt value for testing should be non-negative.\n");

        WntConcentration<2>::Destroy();
    }

    void TestWntInitialisationSetup()
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<WntCellCycleModel*> models;

        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        boost::shared_ptr<AbstractCellProperty> p_stem_type(new StemCellProliferativeType);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            WntCellCycleModel* p_model = new WntCellCycleModel();
            p_model->SetDimension(2);
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            double birth_time = 0.0 - i;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
            models.push_back(p_model);
        }

        // Create the crypt
        MeshBasedCellPopulation<2> crypt(mesh, cells);

        WntConcentration<2>::Instance()->SetType(LINEAR);
        double crypt_length = 1.0;
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);

        // As there is no cell-based simulation we must explicitly initialise the cells
        crypt.InitialiseCells();

        for (AbstractCellPopulation<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            WntCellCycleModel* p_model = static_cast<WntCellCycleModel*>(cell_iter->GetCellCycleModel());
            std::vector<double> proteins = p_model->GetProteinConcentrations();

            if (crypt.GetLocationOfCellCentre(*cell_iter)[1] == 0.0)
            {
                TS_ASSERT_DELTA(proteins[5], 4.975124378109454e-03, 1e-3);
                TS_ASSERT_DELTA(proteins[6]+proteins[7], 6.002649406788524e-01, 1e-3);
                TS_ASSERT_DELTA(proteins[8], 1.00, 1e-3);
            }
            else
            {
                TS_ASSERT_DELTA(crypt.GetLocationOfCellCentre(*cell_iter)[1], 1.0, 1e-12);
                TS_ASSERT_DELTA(proteins[5], 1.000, 1e-3);
                TS_ASSERT_DELTA(proteins[6]+proteins[7], 0.0074, 1e-3);
                TS_ASSERT_DELTA(proteins[8], 0.00, 1e-3);
            }
        }

        // Coverage
        WntConcentration<2>::Instance()->SetConstantWntValueForTesting(5.0);
        c_vector<double, 2> gradient = WntConcentration<2>::Instance()->GetWntGradient(*(crypt.Begin()));
        TS_ASSERT_DELTA(gradient[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(gradient[1], 0.0, 1e-6);

        WntConcentration<2>::Destroy();
    }

    void TestCryptProjectionParameterAAndBGettersAndSetters()
    {
        WntConcentration<2>* p_wnt1 = WntConcentration<2>::Instance();

        p_wnt1->SetCryptProjectionParameterA(0.8);
        p_wnt1->SetCryptProjectionParameterB(1.3);

        WntConcentration<2>* p_wnt2 = WntConcentration<2>::Instance();

        TS_ASSERT_DELTA(p_wnt2->GetCryptProjectionParameterA(), 0.8, 1e-12);
        TS_ASSERT_DELTA(p_wnt2->GetCryptProjectionParameterB(), 1.3, 1e-12);
    }
};

#endif /*TESTWNTCONCENTRATION_HPP_*/
