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


#ifndef TEST1DFIBREWITHHETEROGENEITY_HPP_
#define TEST1DFIBREWITHHETEROGENEITY_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>

#include "TetrahedralMesh.hpp"
#include "Hdf5DataReader.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "MonodomainProblem.hpp"
#include "SimpleStimulus.hpp"

#include "FaberRudy2000.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "NumericFileComparison.hpp"

class HeterogeneousCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    // define a new stimulus
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    HeterogeneousCellFactory()
        : AbstractCardiacCellFactory<1>(),
          mpStimulus(new SimpleStimulus(-600, 0.5))
    {
    }

    HeterogeneousCellFactory(double stimulusMagnitude)
        : AbstractCardiacCellFactory<1>(),
          mpStimulus(new SimpleStimulus(stimulusMagnitude, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        CellFaberRudy2000FromCellML *cell;

        if (this->GetMesh()->GetNode(node)->GetPoint()[0] == 0.0)
        {
            cell = new CellFaberRudy2000FromCellML(this->mpSolver,
                                             mpStimulus);

        }
        else
        {
            cell = new CellFaberRudy2000FromCellML(this->mpSolver,
                                             this->mpZeroStimulus);
        }

        if (this->GetMesh()->GetNode(node)->GetPoint()[0] < 0.3333)
        {
            cell->SetParameter("ScaleFactorGks",0.462);
            cell->SetParameter("ScaleFactorIto",0.0);
        }
        else if (this->GetMesh()->GetNode(node)->GetPoint()[0] < 0.6666)
        {
            cell->SetParameter("ScaleFactorGks",1.154);
            cell->SetParameter("ScaleFactorIto",0.85);
        }
        else //this->mpMesh->GetNode(node)->GetPoint()[0] < 1
        {
            cell->SetParameter("ScaleFactorGks",1.154);
            cell->SetParameter("ScaleFactorIto",1.0);
        }

        return cell;
    }
};


class Test1dFibreWithHeterogeneity : public CxxTest::TestSuite
{
public:
    // Solve on a 1D string of cells, 1cm long with a space step of 0.1mm and heterogeneous cell types.
    void TestFibreHeterogeneity()
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetPdeTimeStep(0.01);
        HeartConfig::Instance()->SetPrintingTimeStep(0.1);
        HeartConfig::Instance()->SetSimulationDuration(300.0);
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_100_elements");
        HeartConfig::Instance()->SetOutputDirectory("FibreWithHeterogeneity");
        HeartConfig::Instance()->SetOutputFilenamePrefix("Monodomain1dWithHeterogeneity");

        HeterogeneousCellFactory cell_factory;
        MonodomainProblem<1> monodomain_problem(&cell_factory);

        monodomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        monodomain_problem.Solve();


        // write out results for node 20 (and 50 and 80)

        Hdf5DataReader results_reader=monodomain_problem.GetDataReader();

        unsigned relevant_nodes[3]={20,50,80};

        for (unsigned i=0; i<3; i++)
        {
            std::vector<double> transmembrane_potential=results_reader.GetVariableOverTime("V", relevant_nodes[i]);
            std::vector<double> time_series = results_reader.GetUnlimitedDimensionValues();

            // Write out the time series for the node at third quadrant
            OutputFileHandler results_handler("FibreWithHeterogeneity", false);
            if (PetscTools::AmMaster())
            {
                OutputFileHandler plot_file_handler("HeterogeneityPlots", false);
                std::stringstream plot_file_name_stream;
                plot_file_name_stream<< "Node_" << relevant_nodes[i] << ".csv";
                out_stream p_plot_file = plot_file_handler.OpenOutputFile(plot_file_name_stream.str());
                for (unsigned data_point = 0; data_point<time_series.size(); data_point++)
                {
                    (*p_plot_file) << time_series[data_point] << "\t" << transmembrane_potential[data_point] << "\n";
                }
                p_plot_file->close();

                std::string file1=plot_file_handler.GetChasteTestOutputDirectory() + "HeterogeneityPlots/" + plot_file_name_stream.str() ;
                std::stringstream file2;
                file2 <<"heart/test/data/HeterogeneityPlots/Node_"<< relevant_nodes[i] << ".csv";
                NumericFileComparison comp(file1, file2.str());
                TS_ASSERT(comp.CompareFiles());
            }
        }

    }
};

#endif /*TEST1DFIBREWITHHETEROGENEITY_HPP_*/
