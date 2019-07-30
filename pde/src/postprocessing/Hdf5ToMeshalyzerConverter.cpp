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

#include "Hdf5ToMeshalyzerConverter.hpp"
#include "MeshalyzerMeshWriter.hpp"
#include "GenericMeshReader.hpp"
#include "UblasCustomFunctions.hpp"
#include "PetscTools.hpp"
#include "Exception.hpp"
#include "ReplicatableVector.hpp"
#include "DistributedVector.hpp"
#include "DistributedVectorFactory.hpp"
#include "Version.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void Hdf5ToMeshalyzerConverter<ELEMENT_DIM,SPACE_DIM>::Write(std::string type)
{

    std::string filename = "";
    if (this->mDatasetNames[this->mOpenDatasetIndex] == "Data")
    {
        filename += this->mFileBaseName + "_";
    }
    filename += type + ".dat";

    out_stream p_file = out_stream(nullptr);
    if (PetscTools::AmMaster())
    {
        p_file = this->mpOutputFileHandler->OpenOutputFile(filename);

        // Check how many digits are to be output in the solution (0 goes to default value of digits)
        if (this->mPrecision != 0)
        {
           p_file->precision(this->mPrecision);
        }
    }

    unsigned num_nodes = this->mpReader->GetNumberOfRows();
    unsigned num_timesteps = this->mpReader->GetUnlimitedDimensionValues().size();

    DistributedVectorFactory factory(num_nodes);

    Vec data = factory.CreateVec();
    ReplicatableVector repl_data(num_nodes);
    for (unsigned time_step=0; time_step<num_timesteps; time_step++)
    {
        this->mpReader->GetVariableOverNodes(data, type, time_step);
        repl_data.ReplicatePetscVector(data);

        assert(repl_data.GetSize() == num_nodes);

        if (PetscTools::AmMaster())
        {
            for (unsigned i=0; i<num_nodes; i++)
            {
                *p_file << repl_data[i] << "\n";
            }
        }
    }
    PetscTools::Destroy(data);
    if (PetscTools::AmMaster())
    {
        std::string comment = "# " + ChasteBuildInfo::GetProvenanceString();
        *p_file << comment;
        p_file->close();
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Hdf5ToMeshalyzerConverter<ELEMENT_DIM,SPACE_DIM>::Hdf5ToMeshalyzerConverter(const FileFinder& rInputDirectory,
                                                                            const std::string& rFileBaseName,
                                                                            AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                                                            bool usingOriginalNodeOrdering,
                                                                            unsigned precision)
    : AbstractHdf5Converter<ELEMENT_DIM,SPACE_DIM>(rInputDirectory, rFileBaseName, pMesh, "output", precision)
{
    do
    {
        std::vector<std::string> variable_names = this->mpReader->GetVariableNames();
        for (unsigned i=0; i<variable_names.size(); i++)
        {
            Write(variable_names[i]);
        }
    }
    while ( this->MoveOntoNextDataset() );

    // Now we might call this class more than once, so we don't always need to write the mesh out.
    // so check to see if it is there already.
    FileFinder test_output("",RelativeTo::ChasteTestOutput);
    std::string output_directory = rInputDirectory.GetRelativePath(test_output) + "/" + this->mRelativeSubdirectory;
    FileFinder mesh_file(output_directory + "/" + rFileBaseName + "_mesh.pts", RelativeTo::ChasteTestOutput);

    if (!mesh_file.IsFile())
    {
        // Write mesh in a suitable form for meshalyzer
        MeshalyzerMeshWriter<ELEMENT_DIM,SPACE_DIM> mesh_writer(output_directory, rFileBaseName + "_mesh", false);

        // Normal case is that the in-memory mesh is converted
        if (!usingOriginalNodeOrdering || !this->mpMesh->IsMeshOnDisk())
        {
            // The second argument tells the writer to not follow original element ordering for performance reasons.
            mesh_writer.WriteFilesUsingMesh(*(this->mpMesh), false);
        }
        else
        {
            // In this case we expect the mesh to have been read in from file
            ///\todo What if the mesh has been scaled, translated or rotated?
            // Note that the next line will throw if the mesh has not been read from file
            std::string original_file = this->mpMesh->GetMeshFileBaseName();
            std::shared_ptr<AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> > p_original_mesh_reader
                = GenericMeshReader<ELEMENT_DIM, SPACE_DIM>(original_file);
            mesh_writer.WriteFilesUsingMeshReader(*p_original_mesh_reader);
        }
    }
    PetscTools::Barrier("Hdf5ToMeshalyzerConverter");
}

// Explicit instantiation
template class Hdf5ToMeshalyzerConverter<1,1>;
template class Hdf5ToMeshalyzerConverter<1,2>;
template class Hdf5ToMeshalyzerConverter<2,2>;
template class Hdf5ToMeshalyzerConverter<1,3>;
template class Hdf5ToMeshalyzerConverter<2,3>;
template class Hdf5ToMeshalyzerConverter<3,3>;
