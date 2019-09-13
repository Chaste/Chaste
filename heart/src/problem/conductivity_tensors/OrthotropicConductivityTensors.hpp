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
#ifndef ORTHOTROPICCONDUCTIVITYTENSORS_HPP_
#define ORTHOTROPICCONDUCTIVITYTENSORS_HPP_

#include <vector>
#include "UblasIncludes.hpp"
#include "AbstractConductivityTensors.hpp"


/**
 *
 * This class provides an abstraction for the definition of constant/non-constant diffusion tensors
 * associated to the different elements of the mesh.
 *
 * After instantiating the class any of SetFibreOrientationFile() or SetNonConstantConductivities()
 * (or both) can be called to implement fibre orientation or heterogeneous conductivity into the
 * tensors, respectively. If none of them is called a constant tensor (with constant conductivities)
 * will be generated for all the elements of the mesh.
 *
 * Init() should be called to actually create the tensors.
 *
 * Initial values for conductivity from "Laminar Arrangement of Ventricular Myocytes Influences Electrical
 * Behavior of the Heart", Hooks et al. 2007
 *
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class OrthotropicConductivityTensors : public AbstractConductivityTensors<ELEMENT_DIM, SPACE_DIM>
{
public:

    /**
     *  Computes the tensors based in all the info set
     *  @param pMesh a pointer to the mesh on which these tensors are to be used
     */
    void Init(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *pMesh);
};

#endif /*ORTHOTROPICCONDUCTIVITYTENSORS_HPP_*/
