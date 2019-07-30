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

#include "AbstractConductivityTensors.hpp"
#include "Exception.hpp"
#include <sstream>

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractConductivityTensors<ELEMENT_DIM,SPACE_DIM>::AbstractConductivityTensors()
    : mpMesh(NULL),
      mUseNonConstantConductivities(false),
      mUseFibreOrientation(false),
      mInitialised(false)
{
    double init_data[]={DBL_MAX, DBL_MAX, DBL_MAX};

    for (unsigned dim=0; dim<SPACE_DIM; dim++)
    {
         mConstantConductivities[dim] = init_data[dim];
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractConductivityTensors<ELEMENT_DIM,SPACE_DIM>::~AbstractConductivityTensors()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractConductivityTensors<ELEMENT_DIM,SPACE_DIM>::SetFibreOrientationFile(const FileFinder &rFibreOrientationFile)
{
    mUseFibreOrientation = true;
    mFibreOrientationFile = rFibreOrientationFile;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractConductivityTensors<ELEMENT_DIM,SPACE_DIM>::SetConstantConductivities(c_vector<double, 1> constantConductivities)
{
    if (SPACE_DIM != 1)
    {
        EXCEPTION("Wrong number of conductivities provided");
    }

    mUseNonConstantConductivities = false;
    mConstantConductivities = constantConductivities;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractConductivityTensors<ELEMENT_DIM,SPACE_DIM>::SetConstantConductivities(c_vector<double, 2> constantConductivities)
{
    if (SPACE_DIM != 2)
    {
        EXCEPTION("Wrong number of conductivities provided");
    }

    mUseNonConstantConductivities = false;
    mConstantConductivities = constantConductivities;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractConductivityTensors<ELEMENT_DIM,SPACE_DIM>::SetConstantConductivities(c_vector<double, 3> constantConductivities)
{
    if (SPACE_DIM != 3)
    {
        EXCEPTION("Wrong number of conductivities provided");
    }

    mUseNonConstantConductivities = false;
    mConstantConductivities = constantConductivities;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractConductivityTensors<ELEMENT_DIM,SPACE_DIM>::SetNonConstantConductivities(std::vector<c_vector<double, SPACE_DIM> >* pNonConstantConductivities)
{
    mUseNonConstantConductivities = true;
    mpNonConstantConductivities = pNonConstantConductivities;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_matrix<double,SPACE_DIM,SPACE_DIM>& AbstractConductivityTensors<ELEMENT_DIM,SPACE_DIM>::operator[](const unsigned global_index)
{
    assert(mInitialised);
    if (global_index >= this->mpMesh->GetNumElements() )
    {
        EXCEPTION("Conductivity tensor requested for element with global_index=" << global_index << ", but there are only " << this->mpMesh->GetNumElements() << " elements in the mesh.");
    }

    if (!mUseNonConstantConductivities && !mUseFibreOrientation)
    {
        return mTensors[0];
    }
    else
    {
        unsigned local_index = mpMesh->SolveElementMapping(global_index); //This will throw if we don't own the element
        return mTensors[local_index];
    }
}

// Explicit instantiation
template class AbstractConductivityTensors<1,1>;
template class AbstractConductivityTensors<1,2>;
template class AbstractConductivityTensors<1,3>;
template class AbstractConductivityTensors<2,2>;
template class AbstractConductivityTensors<2,3>;
template class AbstractConductivityTensors<3,3>;
