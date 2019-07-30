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

#include "VoltageInterpolaterOntoMechanicsMesh.hpp"
#include "Hdf5ToCmguiConverter.hpp"
#include "FineCoarseMeshPair.hpp"
#include "ReplicatableVector.hpp"
#include "HeartConfig.hpp"
#include "Hdf5DataReader.hpp"
#include "PetscTools.hpp"
#include "Hdf5DataWriter.hpp"

template<unsigned DIM>
VoltageInterpolaterOntoMechanicsMesh<DIM>::VoltageInterpolaterOntoMechanicsMesh(
                                     TetrahedralMesh<DIM,DIM>& rElectricsMesh,
                                     QuadraticMesh<DIM>& rMechanicsMesh,
                                     std::vector<std::string>& rVariableNames,
                                     std::string directory,
                                     std::string inputFileNamePrefix)
{
    // Read the data from the HDF5 file
    Hdf5DataReader reader(directory,inputFileNamePrefix);

    unsigned num_timesteps = reader.GetUnlimitedDimensionValues().size();

    // set up the elements and weights for the coarse nodes in the fine mesh
    FineCoarseMeshPair<DIM> mesh_pair(rElectricsMesh, rMechanicsMesh);
    mesh_pair.SetUpBoxesOnFineMesh();
    mesh_pair.ComputeFineElementsAndWeightsForCoarseNodes(true);
    assert(mesh_pair.rGetElementsAndWeights().size()==rMechanicsMesh.GetNumNodes());

    // create and setup a writer
    Hdf5DataWriter* p_writer = new Hdf5DataWriter(*rMechanicsMesh.GetDistributedVectorFactory(),
                                                  directory,
                                                  "voltage_mechanics_mesh",
                                                  false, //don't clean
                                                  false);

    std::vector<int> columns_id;
    for (unsigned var_index = 0; var_index < rVariableNames.size(); var_index++)
    {
        std::string var_name = rVariableNames[var_index];
        columns_id.push_back( p_writer->DefineVariable(var_name,"mV") );
    }

    p_writer->DefineUnlimitedDimension("Time","msecs", num_timesteps);
    p_writer->DefineFixedDimension( rMechanicsMesh.GetNumNodes() );
    p_writer->EndDefineMode();

    assert(columns_id.size() == rVariableNames.size());

    // set up a vector to read into
    DistributedVectorFactory factory(rElectricsMesh.GetNumNodes());
    Vec voltage = factory.CreateVec();
    std::vector<double> interpolated_voltages(rMechanicsMesh.GetNumNodes());
    Vec voltage_coarse = NULL;

    for (unsigned time_step=0; time_step<num_timesteps; time_step++)
    {
        for (unsigned var_index = 0; var_index < rVariableNames.size(); var_index++)
        {
            std::string var_name = rVariableNames[var_index];
            // read
            reader.GetVariableOverNodes(voltage, var_name, time_step);
            ReplicatableVector voltage_repl(voltage);

            // interpolate
            for (unsigned i=0; i<mesh_pair.rGetElementsAndWeights().size(); i++)
            {
                double interpolated_voltage = 0;

                Element<DIM,DIM>& element = *(rElectricsMesh.GetElement(mesh_pair.rGetElementsAndWeights()[i].ElementNum));
                for (unsigned node_index = 0; node_index<element.GetNumNodes(); node_index++)
                {
                    unsigned global_node_index = element.GetNodeGlobalIndex(node_index);
                    interpolated_voltage += voltage_repl[global_node_index]*mesh_pair.rGetElementsAndWeights()[i].Weights(node_index);
                }

                interpolated_voltages[i] = interpolated_voltage;
            }

            if (voltage_coarse != NULL)
            {
                PetscTools::Destroy(voltage_coarse);
            }
            voltage_coarse = PetscTools::CreateVec(interpolated_voltages);
            // write
            p_writer->PutVector(columns_id[var_index], voltage_coarse);
        }
        p_writer->PutUnlimitedVariable(time_step);
        p_writer->AdvanceAlongUnlimitedDimension();
    }

    if (voltage_coarse != NULL)
    {
        PetscTools::Destroy(voltage);
        PetscTools::Destroy(voltage_coarse);
    }

    // delete to flush
    delete p_writer;

    // Convert the new data to CMGUI format.
    // alter the directory in HeartConfig as that is where Hdf5ToCmguiConverter decides
    // where to output
    std::string config_directory = HeartConfig::Instance()->GetOutputDirectory();
    HeartConfig::Instance()->SetOutputDirectory(directory);
    Hdf5ToCmguiConverter<DIM,DIM> converter(FileFinder(directory, RelativeTo::ChasteTestOutput),
                                            "voltage_mechanics_mesh",
                                            &rMechanicsMesh,
                                            false);
    HeartConfig::Instance()->SetOutputDirectory(config_directory);
}

// Explicit instantiation
template class VoltageInterpolaterOntoMechanicsMesh<1>;
template class VoltageInterpolaterOntoMechanicsMesh<2>;
template class VoltageInterpolaterOntoMechanicsMesh<3>;
