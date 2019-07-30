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

#include "UblasCustomFunctions.hpp"
#include "Hdf5ToTxtConverter.hpp"
#include "PetscTools.hpp"
#include "Exception.hpp"
#include "ReplicatableVector.hpp"
#include "DistributedVector.hpp"
#include "DistributedVectorFactory.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "Warnings.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Hdf5ToTxtConverter<ELEMENT_DIM, SPACE_DIM>::Hdf5ToTxtConverter(const FileFinder& rInputDirectory,
                                                               const std::string& rFileBaseName,
                                                               AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh)
    : AbstractHdf5Converter<ELEMENT_DIM,SPACE_DIM>(rInputDirectory, rFileBaseName, pMesh, "txt_output", 0u)
{
    // Make sure that we are never trying to write from an incomplete data HDF5 file
    assert(this->mpReader->GetNumberOfRows() == pMesh->GetNumNodes());

    FileFinder output_directory(this->mRelativeSubdirectory,rInputDirectory);
    OutputFileHandler handler(output_directory);

    unsigned num_nodes = pMesh->GetNumNodes();
    unsigned num_timesteps = this->mpReader->GetUnlimitedDimensionValues().size();

    DistributedVectorFactory* p_factory = pMesh->GetDistributedVectorFactory();
    Vec data = p_factory->CreateVec();
    ReplicatableVector repl_data(num_nodes);

    // Loop over time steps
    for (unsigned time_step=0; time_step<num_timesteps; time_step++)
    {
        // Loop over variables
        for (unsigned var_index=0; var_index<this->mNumVariables; var_index++)
        {
            std::string variable_name = this->mpReader->GetVariableNames()[var_index];

            // Create a .txt file for this time step and this variable
            std::stringstream file_name;
            file_name << rFileBaseName << "_" << variable_name << "_" << time_step << ".txt";
            out_stream p_file = handler.OpenOutputFile(file_name.str());

            this->mpReader->GetVariableOverNodes(data, variable_name, time_step);
            repl_data.ReplicatePetscVector(data);

            assert(repl_data.GetSize() == num_nodes);

            if (PetscTools::AmMaster())
            {
                for (unsigned i=0; i<num_nodes; i++)
                {
                    *p_file << repl_data[i] << "\n";
                }
            }
            p_file->close();
        }
    }

    // Tidy up
    PetscTools::Destroy(data);
}

// Explicit instantiation
template class Hdf5ToTxtConverter<1,1>;
template class Hdf5ToTxtConverter<1,2>;
template class Hdf5ToTxtConverter<2,2>;
template class Hdf5ToTxtConverter<1,3>;
template class Hdf5ToTxtConverter<2,3>;
template class Hdf5ToTxtConverter<3,3>;
