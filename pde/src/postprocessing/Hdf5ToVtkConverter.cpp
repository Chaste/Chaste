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

#include "UblasCustomFunctions.hpp"
#include "Hdf5ToVtkConverter.hpp"
#include "PetscTools.hpp"
#include "Exception.hpp"
#include "ReplicatableVector.hpp"
#include "DistributedVector.hpp"
#include "DistributedVectorFactory.hpp"
#include "VtkMeshWriter.hpp"
#include "GenericMeshReader.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "Warnings.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Hdf5ToVtkConverter<ELEMENT_DIM, SPACE_DIM>::Hdf5ToVtkConverter(std::string inputDirectory,
                                                               std::string fileBaseName,
                                                               AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh,
                                                               bool parallelVtk,
                                                               bool usingOriginalNodeOrdering)
    : AbstractHdf5Converter<ELEMENT_DIM,SPACE_DIM>(inputDirectory, fileBaseName, pMesh, "vtk_output")
{
#ifdef CHASTE_VTK // Requires "sudo aptitude install libvtk5-dev" or similar

    // Write mesh in a suitable form for VTK
    std::string output_directory = inputDirectory + "/" + this->mRelativeSubdirectory;
    VtkMeshWriter<ELEMENT_DIM,SPACE_DIM> vtk_writer(output_directory, fileBaseName, false);

    DistributedVectorFactory* p_factory = pMesh->GetDistributedVectorFactory();

    // Make sure that we are never trying to write from an incomplete data HDF5 file
    assert(this->mpReader->GetNumberOfRows() == pMesh->GetNumNodes());

    DistributedTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* p_distributed_mesh = dynamic_cast<DistributedTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>*>(pMesh);

    unsigned num_nodes = pMesh->GetNumNodes();
    if (parallelVtk)
    {
        // If it's not a distributed mesh, then we might want to give a warning and back-off
        if (p_distributed_mesh == NULL)
        {
            WARNING("Can only write parallel VTK from a DistributedTetrahedralMesh - writing sequential VTK instead");
            parallelVtk = false;
        }

        // If the node ordering flag is set, then we can't do this
        if (usingOriginalNodeOrdering)
        {
            WARNING("Can't write parallel VTK (pvtu) files with original ordering - writing sequential VTK instead");
            parallelVtk = false;
        }

        // Are we now committed to writing .pvtu?
        if (parallelVtk)
        {
           vtk_writer.SetParallelFiles(*pMesh);
           num_nodes = p_distributed_mesh->GetNumLocalNodes();
        }
    }

    Vec data = p_factory->CreateVec();

    unsigned num_timesteps = this->mpReader->GetUnlimitedDimensionValues().size();

    // Loop over time steps
    for (unsigned time_step=0; time_step<num_timesteps; time_step++)
    {
        // Loop over variables
        for (unsigned variable=0; variable<this->mNumVariables; variable++)
        {
            std::string variable_name = this->mpReader->GetVariableNames()[variable];

            // Gets variable at this time step from HDF5 archive
            this->mpReader->GetVariableOverNodes(data, variable_name, time_step);

            std::vector<double> data_for_vtk;
            data_for_vtk.resize(num_nodes);
            std::ostringstream variable_point_data_name;
            variable_point_data_name << variable_name << "_" << std::setw(6) << std::setfill('0') << time_step;

            if (parallelVtk)
            {
                // Parallel VTU files
                double *p_data;
                VecGetArray(data, &p_data);
                for (unsigned index=0; index<num_nodes; index++)
                {
                    data_for_vtk[index]  = p_data[index];
                }
                VecRestoreArray(data, &p_data);
            }
            else
            {
                // One VTU file
                ReplicatableVector repl_data(data);
                for (unsigned index=0; index<num_nodes; index++)
                {
                    data_for_vtk[index] = repl_data[index];
                }
            }
            // Add this variable into the node "point" data
            vtk_writer.AddPointData(variable_point_data_name.str(), data_for_vtk);
        }
    }

    // Tidy up
    PetscTools::Destroy(data);

    // Normally the in-memory mesh is converted
    if (!usingOriginalNodeOrdering)
    {
        vtk_writer.WriteFilesUsingMesh(*(this->mpMesh));
    }
    else
    {
        // In this case we expect the mesh to have been read in from file
        ///\todo What if the mesh has been scaled, translated or rotated?
        // Note that the next line will throw if the mesh has not been read from file
        std::string original_file = this->mpMesh->GetMeshFileBaseName();
        std::auto_ptr<AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> > p_original_mesh_reader
            = GenericMeshReader<ELEMENT_DIM, SPACE_DIM>(original_file);
        vtk_writer.WriteFilesUsingMeshReader(*p_original_mesh_reader);
    }
#endif //CHASTE_VTK
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class Hdf5ToVtkConverter<1,1>;
template class Hdf5ToVtkConverter<1,2>;
template class Hdf5ToVtkConverter<2,2>;
template class Hdf5ToVtkConverter<1,3>;
template class Hdf5ToVtkConverter<2,3>;
template class Hdf5ToVtkConverter<3,3>;
