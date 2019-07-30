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

#include "FibreWriter.hpp"
#include "Version.hpp"

template<unsigned DIM>
FibreWriter<DIM>::FibreWriter(const std::string& rDirectory,
                              const std::string& rBaseName,
                              const bool clearOutputDir)
    : mBaseName(rBaseName),
      mFileIsBinary(false)
{
    mpOutputFileHandler = new OutputFileHandler(rDirectory, clearOutputDir);
}

template<unsigned DIM>
FibreWriter<DIM>::~FibreWriter()
{
    delete mpOutputFileHandler;
}

template<unsigned DIM>
void FibreWriter<DIM>::WriteAllAxi(const std::vector< c_vector<double, DIM> >& fibres)
{
    // Write axi file
    out_stream p_axi_file = OpenFileAndWriteHeader(this->mBaseName + ".axi", fibres.size());

    //Now give the fibre directions
    for (unsigned i=0; i<fibres.size();i++ )
    {
        if (this->mFileIsBinary)
        {
            p_axi_file->write((char*)&fibres[i][0], DIM*sizeof(double));
        }
        else
        {
            for (unsigned j=0; j<DIM; j++)
            {
                *p_axi_file << fibres[i][j] << "\t";
            }
            *p_axi_file <<"\n";
        }
    }
    *p_axi_file << "#\n# " <<  ChasteBuildInfo::GetProvenanceString();
    p_axi_file->close();
}

template<unsigned DIM>
void FibreWriter<DIM>::WriteAllOrtho(const std::vector< c_vector<double, DIM> >& fibres,
                                     const std::vector< c_vector<double, DIM> >& second,
                                     const std::vector< c_vector<double, DIM> >& third)
{
    assert(fibres.size() == second.size());
    assert(second.size() == third.size());
    // Write ortho file
    out_stream p_file = OpenFileAndWriteHeader(this->mBaseName + ".ortho", fibres.size());

    //Now give the fibre directions
    for (unsigned i=0; i<fibres.size();i++ )
    {
        if (this->mFileIsBinary)
        {
            //The binary file is row-major
            p_file->write((char*)&fibres[i][0], DIM*sizeof(double));
            p_file->write((char*)&second[i][0], DIM*sizeof(double));
            p_file->write((char*)&third[i][0], DIM*sizeof(double));
        }
        else
        {
            //The ascii file is row-major
            for (unsigned j=0; j<DIM; j++)
            {
                *p_file << fibres[i][j] << "\t";
            }
            for (unsigned j=0; j<DIM; j++)
            {
                *p_file << second[i][j] << "\t";
            }
            for (unsigned j=0; j<DIM; j++)
            {
                *p_file << third[i][j] << "\t";
            }
            *p_file <<"\n";
        }
    }
    *p_file << "#\n# " <<  ChasteBuildInfo::GetProvenanceString();
    p_file->close();
}

template<unsigned DIM>
out_stream FibreWriter<DIM>::OpenFileAndWriteHeader(const std::string& rFileName, unsigned numItems)
{
    out_stream p_fibre_file = this->mpOutputFileHandler->OpenOutputFile(rFileName);

    //Header line says how many fibres are coming...
    *p_fibre_file << numItems;
    //... and whether the fibres are binary
    if (this->mFileIsBinary)
    {
        *p_fibre_file << "\tBIN\n";
    }
    else
    {
        *p_fibre_file << "\n";
    }
    return p_fibre_file;
}

template<unsigned DIM>
void FibreWriter<DIM>::SetWriteFileAsBinary()
{
    mFileIsBinary = true;
}

// Explicit instantiation
template class FibreWriter<1>;
template class FibreWriter<2>;
template class FibreWriter<3>;
