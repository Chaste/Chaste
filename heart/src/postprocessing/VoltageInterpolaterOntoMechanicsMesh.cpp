/*

Copyright (C) University of Oxford, 2005-2011

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

#include "VoltageInterpolaterOntoMechanicsMesh.hpp"

template<unsigned DIM>
VoltageInterpolaterOntoMechanicsMesh<DIM>::VoltageInterpolaterOntoMechanicsMesh(
                                     TetrahedralMesh<DIM,DIM>& rElectricsMesh,
                                     QuadraticMesh<DIM>& rMechanicsMesh,
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

    p_writer->DefineFixedDimension( rMechanicsMesh.GetNumNodes() );
    int voltage_column_id = p_writer->DefineVariable("V","mV");
    p_writer->DefineUnlimitedDimension("Time","msecs");
    p_writer->EndDefineMode();

    // set up a vector to read into
    DistributedVectorFactory factory(rElectricsMesh.GetNumNodes());
    Vec voltage = factory.CreateVec();
    std::vector<double> interpolated_voltages(rMechanicsMesh.GetNumNodes());
    Vec voltage_coarse = NULL;

    for(unsigned time_step=0; time_step<num_timesteps; time_step++)
    {
        // read
        reader.GetVariableOverNodes(voltage, "V", time_step);
        ReplicatableVector voltage_repl(voltage);

        // interpolate
        for(unsigned i=0; i<mesh_pair.rGetElementsAndWeights().size(); i++)
        {
            double interpolated_voltage = 0;

            Element<DIM,DIM>& element = *(rElectricsMesh.GetElement(mesh_pair.rGetElementsAndWeights()[i].ElementNum));
            for(unsigned node_index = 0; node_index<element.GetNumNodes(); node_index++)
            {
                unsigned global_node_index = element.GetNodeGlobalIndex(node_index);
                interpolated_voltage += voltage_repl[global_node_index]*mesh_pair.rGetElementsAndWeights()[i].Weights(node_index);
            }

            interpolated_voltages[i] = interpolated_voltage;
        }

        if(voltage_coarse!=NULL)
        {
            VecDestroy(voltage_coarse);
        }
        voltage_coarse = PetscTools::CreateVec(interpolated_voltages);

        // write
        p_writer->PutUnlimitedVariable(time_step);
        p_writer->PutVector(voltage_column_id, voltage_coarse);
        p_writer->AdvanceAlongUnlimitedDimension();
    }

    if(voltage_coarse!=NULL)
    {
        VecDestroy(voltage);
        VecDestroy(voltage_coarse);
    }

    // delete to flush
    delete p_writer;

    // Convert the new data to CMGUI format.
    // alter the directory in HeartConfig as that is where Hdf5ToCmguiConverter decides
    // where to output
    std::string config_directory = HeartConfig::Instance()->GetOutputDirectory();
    HeartConfig::Instance()->SetOutputDirectory(directory);
    Hdf5ToCmguiConverter<DIM,DIM> converter(directory,
                                            "voltage_mechanics_mesh",
                                            &rMechanicsMesh,
                                            false);
    HeartConfig::Instance()->SetOutputDirectory(config_directory);
}

///////////////////////////////////////////
// explicit instantiation
///////////////////////////////////////////

template class VoltageInterpolaterOntoMechanicsMesh<1>;
template class VoltageInterpolaterOntoMechanicsMesh<2>;
template class VoltageInterpolaterOntoMechanicsMesh<3>;
