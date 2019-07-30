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


#ifndef ABSTRACTCORRECTIONTERMASSEMBLER_HPP_
#define ABSTRACTCORRECTIONTERMASSEMBLER_HPP_


#include "AbstractCardiacFeVolumeIntegralAssembler.hpp"
#include "AbstractCardiacTissue.hpp"

/**
 * A parent class for MonodomainCorrectionTermAssembler and BidomainCorrectionTermAssembler,
 * used for state variable interpolation (SVI).
 */
template<unsigned ELEM_DIM,unsigned SPACE_DIM,unsigned PROBLEM_DIM>
class AbstractCorrectionTermAssembler
    : public AbstractCardiacFeVolumeIntegralAssembler<ELEM_DIM,SPACE_DIM,PROBLEM_DIM,true,false,CARDIAC>
{
protected:
    /** Ionic current to be interpolated from cache */
    double mIionicInterp;

    /** State variables interpolated onto quadrature point */
    std::vector<double> mStateVariablesAtQuadPoint;

    /**
     * Resets interpolated state variables and ionic current.
     */
    void ResetInterpolatedQuantities( void );

    /**
     *  Vector of bools, one bool per element, representing whether each
     *  element can do SVI. If the element has different cell models at
     *  each node, or has been designated a bath element, then it cannot
     *  do SVI.
     */
    std::vector<bool> mElementsCanDoSvi;

    /**
     * Interpolates state variables and ionic current.
     *
     * @param phiI
     * @param pNode
     */
    void IncrementInterpolatedQuantities(double phiI, const Node<SPACE_DIM>* pNode);

    /**
     * @return true if we should assemble the correction term for this element.
     * Checks if there is a sufficiently steep ionic current gradient to make the expense worthwhile, by checking
     * if the maximum difference between nodal ionic currents is greater than 1 uA/cm^2^.
     *
     * @param rElement  the element to test
     */
    bool ElementAssemblyCriterion(Element<ELEM_DIM,SPACE_DIM>& rElement);

public:

    /**
     * Constructor.
     *
     * @param pMesh  pointer to the mesh
     * @param pTissue  pointer to the cardiac tissue
     */
    AbstractCorrectionTermAssembler(AbstractTetrahedralMesh<ELEM_DIM,SPACE_DIM>* pMesh,
                                    AbstractCardiacTissue<ELEM_DIM,SPACE_DIM>* pTissue);
};

#endif /*ABSTRACTCORRECTIONTERMASSEMBLER_HPP_*/
