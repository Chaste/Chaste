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

    if (!mUseNonConstantConductivities && !mUseFibreOrientation)
    {
        return mTensors[0];
    }
    else
    {
        assert(global_index < this->mpMesh->GetNumElements());
        unsigned local_index = mpMesh->SolveElementMapping(global_index); //This will throw if we don't own the element
        return mTensors[local_index];
    }
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class AbstractConductivityTensors<1,1>;
template class AbstractConductivityTensors<1,2>;
template class AbstractConductivityTensors<1,3>;
template class AbstractConductivityTensors<2,2>;
template class AbstractConductivityTensors<2,3>;
template class AbstractConductivityTensors<3,3>;
