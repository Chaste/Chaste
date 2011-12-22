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


#ifndef ABSTRACTCORRECTIONTERMASSEMBLER_HPP_
#define ABSTRACTCORRECTIONTERMASSEMBLER_HPP_


#include "AbstractCardiacFeVolumeIntegralAssembler.hpp"
#include "AbstractCardiacTissue.hpp"

/**
 * A parent class for MonodomainCorrectionTermAssembler and BidomainCorrectionTermAssembler
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
     *  Vector of bools, one bool per element, saying whether that
     *  element has identical cell models at each node. If this
     *  is not the case, SVI is certainly not posssible in this element
     */
    std::vector<bool> mElementsHasIdenticalCellModels;

    /**
     * Interpolates state variables and ionic current.
     *
     * @param phiI
     * @param pNode
     */
    void IncrementInterpolatedQuantities(double phiI, const Node<SPACE_DIM>* pNode);


    /**
     * Determine whether to assemble the correction term for this element.
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
     * @param numQuadPoints  number of quadrature points
     */
    AbstractCorrectionTermAssembler(AbstractTetrahedralMesh<ELEM_DIM,SPACE_DIM>* pMesh,
                                    AbstractCardiacTissue<ELEM_DIM,SPACE_DIM>* pTissue,
                                    unsigned numQuadPoints = 2);
};



#endif /*ABSTRACTCORRECTIONTERMASSEMBLER_HPP_*/
