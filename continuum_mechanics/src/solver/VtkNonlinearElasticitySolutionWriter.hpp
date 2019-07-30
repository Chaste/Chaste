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

#ifndef VTKNONLINEARELASTICITYSOLUTIONWRITER_HPP_
#define VTKNONLINEARELASTICITYSOLUTIONWRITER_HPP_

#include "AbstractNonlinearElasticitySolver.hpp"

/**
 *  Class for write mechanics solutions to .vtu file (for visualisation in Paraview), including
 *  displacement, pressure if incompressible simulation, different strains, and (in future) stresses.
 */
template<unsigned DIM>
class VtkNonlinearElasticitySolutionWriter
{
friend class TestVtkNonlinearElasticitySolutionWriter;

private:
    /** Pointer to the mechanics solver which performed the calculation */
    AbstractNonlinearElasticitySolver<DIM>* mpSolver;
    /** Whether to write strains for each element */
    bool mWriteElementWiseStrains;
    /** What type of strain to write for each element, from: F = dx/dX, C = F^T F, E = 1/2 (C-I) */
    StrainType mElementWiseStrainType;

    /** Tensor data to be written to the .vtu file. This is a member variable only for testing reasons. */
    std::vector<c_matrix<double,DIM,DIM> > mTensorData;


    //// For future..
    //    bool mWriteNodewiseStresses;
    //    StressType mNodeWiseStressType;

public:

    /**
     *  Constructor
     *  @param rSolver mechanics solver which performed the calculation
     */
    VtkNonlinearElasticitySolutionWriter(AbstractNonlinearElasticitySolver<DIM>& rSolver)
        : mpSolver(&rSolver),
          mWriteElementWiseStrains(false)
    {
    }

    /**
     *  Set write strains for each element. Can write any of: F = dx/dX, C = F^T F, E = 1/2 (C-I)
     *  @param strainType Which strain to write, choose one of: DEFORMATION_GRADIENT_F, DEFORMATION_TENSOR_C, LAGRANGE_STRAIN_E
     */
    void SetWriteElementWiseStrains(StrainType strainType)
    {
        mWriteElementWiseStrains = true;
        mElementWiseStrainType = strainType;
    }

    /** Write the .vtu file */
    void Write();
};

#endif // VTKNONLINEARELASTICITYSOLUTIONWRITER_HPP_
