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

#include "AbstractHdf5Converter.hpp"
#include "Version.hpp"


/*
 * Operator function to be called by H5Literate [HDF5 1.8.x] or H5Giterate [HDF5 1.6.x] (in TestListingDatasetsInAnHdf5File).
 */
herr_t op_func (hid_t loc_id,
                const char *name,
                const H5L_info_t *info,
                void *operator_data);


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractHdf5Converter<ELEMENT_DIM, SPACE_DIM>::AbstractHdf5Converter(const FileFinder& rInputDirectory,
                                                                     const std::string& rFileBaseName,
                                                                     AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh,
                                                                     const std::string& rSubdirectoryName,
                                                                     unsigned precision)
    : mrH5Folder(rInputDirectory),
      mFileBaseName(rFileBaseName),
      mOpenDatasetIndex(UNSIGNED_UNSET),
      mpMesh(pMesh),
      mRelativeSubdirectory(rSubdirectoryName),
      mPrecision(precision)
{
    GenerateListOfDatasets(mrH5Folder, mFileBaseName);

    // Create new directory in which to store everything
    FileFinder sub_directory(mRelativeSubdirectory, mrH5Folder);
    mpOutputFileHandler = new OutputFileHandler(sub_directory, false);

    MoveOntoNextDataset();
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
        if (mDatasetNames[mOpenDatasetIndex]=="Data")
        {
           time_info_filename = mFileBaseName + "_times.info";
        }
        else
        {
           time_info_filename = mDatasetNames[mOpenDatasetIndex] + "_times.info";
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
    delete mpOutputFileHandler;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::string AbstractHdf5Converter<ELEMENT_DIM,SPACE_DIM>::GetSubdirectory()
{
    return mRelativeSubdirectory;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractHdf5Converter<ELEMENT_DIM,SPACE_DIM>::MoveOntoNextDataset()
{
    // If we are already at the end just return false.
    if (mDatasetNames.size() == mOpenDatasetIndex+1u)
    {
        return false;
    }

    // If we haven't read anything yet, start at the beginning, otherwise increment by one.
    if (mOpenDatasetIndex==UNSIGNED_UNSET)
    {
        mOpenDatasetIndex = 0u;
    }
    else
    {
        mOpenDatasetIndex++;
    }

    // Store directory, mesh and filenames and create the reader
    mpReader.reset(new Hdf5DataReader(mrH5Folder, mFileBaseName, mDatasetNames[mOpenDatasetIndex]));

    // Check the data file for basic validity
    std::vector<std::string> variable_names = mpReader->GetVariableNames();
    mNumVariables = variable_names.size();

    if (mpReader->GetNumberOfRows() != mpMesh->GetNumNodes())
    {
        delete mpOutputFileHandler;
        EXCEPTION("Mesh and HDF5 file have a different number of nodes");
    }
    WriteInfoFile();

    return true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractHdf5Converter<ELEMENT_DIM,SPACE_DIM>::GenerateListOfDatasets(const FileFinder& rH5Folder,
                                                                          const std::string& rFileName)
{
    /*
     * Open file.
     */
    std::string file_name = rH5Folder.GetAbsolutePath() + rFileName + ".h5";
    hid_t file = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    /*
     * Begin HDF5 iteration, calls a method that populates mDatasetNames.
     */
    H5Literate(file, H5_INDEX_NAME, H5_ITER_NATIVE, nullptr, op_func, &mDatasetNames);

    H5Fclose(file);

    // Remove datasets that end in "_Unlimited", as these are paired up with other ones!
    std::string ending = "_Unlimited";

    // Strip off the independent variables from the list
    std::vector<std::string>::iterator iter;
    for (iter = mDatasetNames.begin(); iter != mDatasetNames.end(); )
    {
        // If the dataset name is "Time" OR ...
        // it is longer than the ending we are looking for ("_Unlimited") ...
        // ... AND it ends with the string we are looking for,
        // then erase it.
        if ((*(iter) == "Time") ||
            ((iter->length() > ending.length()) &&
            (0 == iter->compare(iter->length() - ending.length(), ending.length(), ending))))
        {
            iter = mDatasetNames.erase(iter);
        }
        else
        {
            ++iter;
        }
    }
}

/*
 * HDF5 Operator function.
 *
 * Puts the name of the objects (in this case 'datasets')
 * in an HDF5 file into a std::vector for us to use for
 * iterating over the file.
 *
 * This was based on a couple of HDF5 example files.
 */
herr_t op_func (hid_t loc_id, const char *name,
                const H5L_info_t *info,
                void *operator_data)
{
    std::vector<std::string>* p_dataset_names = static_cast<std::vector< std::string > * >(operator_data);

    /*
     * Get type of the object and display its name and type.
     * The name of the object is passed to this function by
     * the Library.
     */
    H5O_info_t infobuf;
    H5Oget_info_by_name (loc_id, name, &infobuf, H5P_DEFAULT);
    switch (infobuf.type)
    {
//          case H5O_TYPE_GROUP:
//              printf ("  Group: %s\n", name);
//              break;
        case H5O_TYPE_DATASET:
            p_dataset_names->push_back(name);
            break;
//          case H5O_TYPE_NAMED_DATATYPE:
//              printf ("  Datatype: %s\n", name);
//              break;
        default:
            NEVER_REACHED;
            // If you ever do reach here, it means that an HDF5 file you are trying to convert contains
            // something other than a 'Dataset', which is the usual data structure we write out in Chaste.
            // The above commented out lines should help you figure out what it is, and how it got there.
    }
    return 0;
}

// Explicit instantiation
template class AbstractHdf5Converter<1,1>;
template class AbstractHdf5Converter<1,2>;
template class AbstractHdf5Converter<2,2>;
template class AbstractHdf5Converter<1,3>;
template class AbstractHdf5Converter<2,3>;
template class AbstractHdf5Converter<3,3>;
