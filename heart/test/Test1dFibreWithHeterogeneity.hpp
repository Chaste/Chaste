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
