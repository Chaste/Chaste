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

#include "OrthotropicConductivityTensors.hpp"
#include "Exception.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OrthotropicConductivityTensors<ELEMENT_DIM, SPACE_DIM>::Init(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *pMesh)
{
    this->mpMesh = pMesh;

    if (!this->mUseNonConstantConductivities && !this->mUseFibreOrientation)
    {
        // Constant tensor for every element
        c_matrix<double, SPACE_DIM, SPACE_DIM> conductivity_matrix(zero_matrix<double>(SPACE_DIM,SPACE_DIM));

        for (unsigned dim=0; dim<SPACE_DIM; dim++)
        {
            assert(this->mConstantConductivities(dim) != DBL_MAX);
            conductivity_matrix(dim,dim) = this->mConstantConductivities(dim);
        }

        this->mTensors.push_back(conductivity_matrix);
    }
    else
    {
        c_matrix<double,SPACE_DIM,SPACE_DIM> orientation_matrix((identity_matrix<double>(SPACE_DIM)));

        if (this->mUseFibreOrientation)
        {
            // open file
            this->mFileReader.reset(new FibreReader<SPACE_DIM>(this->mFibreOrientationFile, ORTHO));
            if (this->mFileReader->GetNumLinesOfData() != this->mpMesh->GetNumElements())
            {
                EXCEPTION("The size of the fibre file does not match the number of elements in the mesh");
            }
        }

        if (this->mUseNonConstantConductivities)
        {
            if (this->mpNonConstantConductivities->size() != this->mpMesh->GetNumLocalElements())
            {
                EXCEPTION("The size of the conductivities vector does not match the number of elements in the mesh");
            }
        }

        // reserve() allocates all the memory at once, more efficient than relying
        // on the automatic reallocation scheme.
        this->mTensors.reserve(this->mpMesh->GetNumLocalElements());

        c_matrix<double, SPACE_DIM, SPACE_DIM> conductivity_matrix(zero_matrix<double>(SPACE_DIM,SPACE_DIM));

        unsigned local_element_index = 0;

        for (typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator it = this->mpMesh->GetElementIteratorBegin();
             it != this->mpMesh->GetElementIteratorEnd();
             ++it)
        {
            /*
             *  For every element of the mesh we compute its tensor like (from
             * "Laminar Arrangement of Ventricular Myocytes Influences Electrical
             * Behavior of the Heart, DA Hooks et al. 2007):
             *
             *
             *                         [g_f  0   0 ] [a_f']
             *  tensor = [a_f a_l a_n] [ 0  g_l  0 ] [a_l']
             *                         [ 0   0  g_n] [a_n']
             *
             *              [x_i]
             *  where a_i = [y_i], i={f,l,n} are read from the fibre orientation file and
             *              [z_i]
             *
             *  g_f = fibre/longitudinal conductivity (constant or element specific)
             *  g_l = laminar/transverse conductivity (constant or element specific)
             *  g_n = normal conductivity (constant or element specific)
             *
             */
            if (this->mUseNonConstantConductivities)
            {
                for (unsigned dim=0; dim<SPACE_DIM; dim++)
                {
                    conductivity_matrix(dim,dim) = (*this->mpNonConstantConductivities)[local_element_index][dim];
                }
            }
            else
            {
                for (unsigned dim=0; dim<SPACE_DIM; dim++)
                {
                    assert(this->mConstantConductivities(dim) != DBL_MAX);
                    conductivity_matrix(dim,dim) = this->mConstantConductivities(dim);
                }
            }

            if (this->mUseFibreOrientation)
            {
                unsigned current_fibre_global_index = it->GetIndex();
                this->mFileReader->GetFibreSheetAndNormalMatrix(current_fibre_global_index, orientation_matrix);
            }

            c_matrix<double,SPACE_DIM,SPACE_DIM> temp;
            noalias(temp) = prod(orientation_matrix, conductivity_matrix);
            this->mTensors.push_back( prod(temp, trans(orientation_matrix) ) );

            local_element_index++;
        }
        assert(this->mTensors.size() == this->mpMesh->GetNumLocalElements());
        assert(this->mTensors.size() == local_element_index);

        if (this->mUseFibreOrientation)
        {
            // Close fibre file
            this->mFileReader.reset();
        }
    }

    this->mInitialised = true;
}

// Explicit instantiation
template class OrthotropicConductivityTensors<1,1>;
template class OrthotropicConductivityTensors<1,2>;
template class OrthotropicConductivityTensors<1,3>;
template class OrthotropicConductivityTensors<2,2>;
template class OrthotropicConductivityTensors<2,3>;
template class OrthotropicConductivityTensors<3,3>;
