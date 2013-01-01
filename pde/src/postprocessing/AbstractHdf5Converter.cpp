/*

Copyright (c) 2005-2013, University of Oxford.
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

#include "AbstractHdf5Converter.hpp"
#include "Version.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractHdf5Converter<ELEMENT_DIM, SPACE_DIM>::AbstractHdf5Converter(std::string inputDirectory,
                                                                     std::string fileBaseName,
                                                                     AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh,
                                                                     std::string subdirectoryName,
                                                                     std::string datasetName)
    : mFileBaseName(fileBaseName),
      mDatasetName(datasetName),
      mpMesh(pMesh),
      mRelativeSubdirectory(subdirectoryName)
{
    // Store directory, mesh and filenames and create the reader
    mpReader = new Hdf5DataReader(inputDirectory, mFileBaseName, true, mDatasetName);

    // Create new directory in which to store everything
    mpOutputFileHandler = new OutputFileHandler(inputDirectory + "/" + mRelativeSubdirectory, false);

    // Check the data file for basic validity
    std::vector<std::string> variable_names = mpReader->GetVariableNames();
    mNumVariables = variable_names.size();

    if (mpReader->GetNumberOfRows() != mpMesh->GetNumNodes())
    {
        delete mpReader;
        delete mpOutputFileHandler;
        EXCEPTION("Mesh and HDF5 file have a different number of nodes");
    }

    WriteInfoFile();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractHdf5Converter<ELEMENT_DIM, SPACE_DIM>::WriteInfoFile()
{
    // Note that we don't want the child processes to write info files
    if (PetscTools::AmMaster())
    {
        std::string time_info_filename;

        // If the dataset is just "Data" then we will leave the original filename as it is (to avoid confusion!)
        // If the dataset is a new variant like "Postprocessing" then we will put the dataset name in the output.
        if (mDatasetName=="Data")
        {
           time_info_filename = mFileBaseName + "_times.info";
        }
        else
        {
           time_info_filename = mFileBaseName + "_" + mDatasetName + "_times.info";
        }
        out_stream p_file = mpOutputFileHandler->OpenOutputFile(time_info_filename);

        std::vector<double> time_values = mpReader->GetUnlimitedDimensionValues();
        unsigned num_timesteps = time_values.size();
        double first_timestep = time_values.front();
        double last_timestep = time_values.back();

        double timestep = num_timesteps > 1 ? time_values[1] - time_values[0] : DOUBLE_UNSET;

        *p_file << "Number of timesteps " << num_timesteps << std::endl;
        *p_file << "timestep " << timestep << std::endl;
        *p_file << "First timestep " << first_timestep << std::endl;
        *p_file << "Last timestep " << last_timestep << std::endl;
        *p_file << ChasteBuildInfo::GetProvenanceString();

        p_file->close();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractHdf5Converter<ELEMENT_DIM,SPACE_DIM>::~AbstractHdf5Converter()
{
    delete mpReader;
    delete mpOutputFileHandler;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::string AbstractHdf5Converter<ELEMENT_DIM,SPACE_DIM>::GetSubdirectory()
{
    return mRelativeSubdirectory;
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class AbstractHdf5Converter<1,1>;
template class AbstractHdf5Converter<1,2>;
template class AbstractHdf5Converter<2,2>;
template class AbstractHdf5Converter<1,3>;
template class AbstractHdf5Converter<2,3>;
template class AbstractHdf5Converter<3,3>;
