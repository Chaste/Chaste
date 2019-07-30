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
#ifndef ABSTRACTPERELEMENTWRITER_HPP_
#define ABSTRACTPERELEMENTWRITER_HPP_

#include "AbstractTetrahedralMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "OutputFileHandler.hpp"
#include "Version.hpp"

/**
 * An abstract writer class for writing stuff on a "per element" basis.
 * This class will "visit" all the locally owned elements and concentrate
 * data back to the master.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned DATA_SIZE>
class AbstractPerElementWriter
{
private:
    bool mFileIsBinary;  /**< Whether all data is to be written as binary*/

protected:
    /**
     * The mesh. Set by the constructor.
     */
    AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* mpMesh;

    /**
     * The output file (only valid on master process).
     * Set by the public method WriteData and used by WriteElementOnMaster
     */
    out_stream mpMasterFile;

    /**
     * How to associate an element with some data
     * Must be over-ridden by the derived class.
     *
     * @param pElement  a locally-owned element for which to calculate or lookup some data
     * @param localElementIndex the index of pElement in the local vector.  Used in subclasses which look up data from a separate structure ordered by local indices
     * @param rData  the double-precision data to write to file (output from the method)
     */
    virtual void Visit(Element<ELEMENT_DIM, SPACE_DIM>* pElement,
                       unsigned localElementIndex,
                       c_vector<double, DATA_SIZE>& rData)=0;

    /**
     * How to write an element's worth of data to the file.
     * By default writes tab-separated data to a single line, but can be over-ridden.
     * This is only called by the master process.
     *
     * @param rData  the double-precision data to write to file
     */
    virtual void WriteElementOnMaster(const c_vector<double, DATA_SIZE>& rData)
    {
        if (mFileIsBinary)
        {
            //The binary file is row-major
            mpMasterFile->write((char*)&rData[0], DATA_SIZE*sizeof(double));
        }
        else
        {
            for (unsigned i=0; i<DATA_SIZE; i++)
            {
                (*mpMasterFile) << rData[i] << "\t";
            }
            (*mpMasterFile)<<"\n";
        }
    }

    /**
     * How to write the header information to the file.
     * By default writes nothing.
     * This is only called by the master process.
     * This should NOT end the line (eg: \n or std::endl)
     * as we need to say whether the file is binary or not.
     */
    virtual void WriteHeaderOnMaster()
    {
    }

    /**
     * Method that can be overridden to do any pre-calculations necessary to
     * write out data.
     *
     * @param rOutputDirectory  The folder data is going to be written into (mostly for debugging to be written into).
     */
    virtual void PreWriteCalculations(OutputFileHandler& rOutputDirectory)
    {
    }

public:

    /**
     * Constructor
     *
     * @param pMesh  The mesh whose elements we are going to write out data for.
     */
    AbstractPerElementWriter(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh)
     : mFileIsBinary(false),
       mpMesh(pMesh),
       mpMasterFile(nullptr)
    {

    }

    /**
     * Writes data about each element in parallel
     * Data about each element is retrieved by the Visit() method.
     * Writing is done by the master process using the WriteElement() method.
     * Any element not owned by the master is communicated by the unique designated owner.
     *
     * MUST BE CALLED IN PARALLEL.
     *
     * @param rHandler  specify the directory in which to place the output file
     * @param rFileName  the file name
     */
    void WriteData(OutputFileHandler& rHandler, const std::string& rFileName)
    {
        PreWriteCalculations(rHandler);

        c_vector<double, DATA_SIZE> data;
        if (PetscTools::AmMaster())
        {
            mpMasterFile = rHandler.OpenOutputFile(rFileName);
            MPI_Status status;
            status.MPI_ERROR = MPI_SUCCESS; //For MPICH2
            WriteHeaderOnMaster();
            // say whether the fibres are binary in the header line (not ended yet!)
            if (mFileIsBinary)
            {
                *mpMasterFile << "\tBIN\n";
            }
            else
            {
                *mpMasterFile << "\n";
            }
            // The master process needs to keep track of both the global element list
            // (so that all elements are concentrated) and the local element list (so that
            // a local element index can be applied
            unsigned local_element_index=0u; //Data invariant associates iter with local_element_index
            typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator iter = mpMesh->GetElementIteratorBegin();

            for (unsigned global_element_index=0; global_element_index<mpMesh->GetNumElements(); global_element_index++)
            {
                if (mpMesh->CalculateDesignatedOwnershipOfElement(global_element_index))
                {
                    Element<ELEMENT_DIM,SPACE_DIM>* p_elem = mpMesh->GetElement(global_element_index);
                    //Spool forward in the local vector of elements
                    while (&(*iter) != p_elem)
                    {
                        ++iter;
                        local_element_index++;
                        assert(iter != mpMesh->GetElementIteratorEnd());
                    }
                    //Master owns this process and can write it directly
                    Visit(p_elem, local_element_index, data);
                }
                else
                {
                    //Data must come from a remote process
                    MPI_Recv(&data[0], DATA_SIZE, MPI_DOUBLE, MPI_ANY_SOURCE, global_element_index, PETSC_COMM_WORLD, &status);
                }
                WriteElementOnMaster(data);
            }
            *mpMasterFile << "# "<<ChasteBuildInfo::GetProvenanceString();
            mpMasterFile->close();
        }
        else
        {
            //Not master process
            unsigned previous_index = 0u; //Used to check that the indices are monotone
            unsigned local_index = 0u;
            for (typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator iter = mpMesh->GetElementIteratorBegin();
                          iter != mpMesh->GetElementIteratorEnd();
                          ++iter, local_index++)
            {
                unsigned element_index = iter->GetIndex();
                //Check monotonicity
                if (previous_index>0u)
                {
                    assert(element_index > previous_index);
                }
                previous_index = element_index;
                if (mpMesh->CalculateDesignatedOwnershipOfElement(element_index))
                {
                    // The master needs to know about this one.
                    Visit(&(*iter), local_index, data);
                    /// \todo See if this can be speeded up with #2351.
                    MPI_Ssend(&data[0], DATA_SIZE, MPI_DOUBLE, 0, element_index, PETSC_COMM_WORLD);//Tag with element_index
                }
            }

        }
    }

    /**
     * Switch to write binary fibre file
     *
     * (set to write ascii files in the constructor)
     *
     * @param binary  Whether to write as binary (defaults to true).
     */
     void SetWriteFileAsBinary(bool binary=true)
     {
         mFileIsBinary = binary;
     }

    /**
     * Empty virtual destructor for abstract class
     */
    virtual ~AbstractPerElementWriter()
    {
    }
};


#endif /*ABSTRACTPERELEMENTWRITER_HPP_*/
