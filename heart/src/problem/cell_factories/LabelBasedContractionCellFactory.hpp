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

#ifndef LABELBASEDCONTRACTIONCELLFACTORY_HPP_
#define LABELBASEDCONTRACTIONCELLFACTORY_HPP_

#include "AbstractContractionCellFactory.hpp"
#include "ContractionModelName.hpp"
#include "NonPhysiologicalContractionModel.hpp"
#include "Kerchoffs2003ContractionModel.hpp"
#include "Nash2004ContractionModel.hpp"
#include "NhsModelWithBackwardSolver.hpp"
#include "ConstantActiveTension.hpp"

/**
 * This class provides homogeneous cells of a certain pre-specified type
 * (from the list in ContractionModelName.hpp) across the whole mechanics
 * mesh.
 */
template <unsigned DIM>
class LabelBasedContractionCellFactory : public AbstractContractionCellFactory<DIM>
{
private:
    /** The name of the contraction model to supply */
    ContractionModelName mContractionModelName;

public:
    /**
     * Constructor
     *
     * @param contractionModelName  The name of the model to supply to the whole mesh.
     */
    LabelBasedContractionCellFactory(ContractionModelName contractionModelName)
        : AbstractContractionCellFactory<DIM>(),
          mContractionModelName(contractionModelName)
    {
    }

    /**
     * Main method that is called by the solver to assign a model to each quadrature point.
     *
     * @param pElement  Pointer to the element in the mechanics mesh
     * @return  A contraction model for this element.
     */
    AbstractContractionModel* CreateContractionCellForElement(Element<DIM, DIM>* pElement)
    {
        switch(mContractionModelName)
        {
            case CONSTANT:
            {
                return new ConstantActiveTension;
                break;
            }
            case NONPHYSIOL1:
            case NONPHYSIOL2:
            case NONPHYSIOL3:
            {
                unsigned option = (mContractionModelName==NONPHYSIOL1 ? 1 : (mContractionModelName==NONPHYSIOL2? 2 : 3));
                return new NonPhysiologicalContractionModel(option);
                break;
            }
            case NASH2004: //stretch dependent, will this work with explicit??
            {
                return new Nash2004ContractionModel;
                break;
            }
            case NHS:
            {
                return new NhsModelWithBackwardSolver;
                break;
            }
            case KERCHOFFS2003: //stretch dependent, will this work with explicit? Answer: can be unstable
            {
                return new Kerchoffs2003ContractionModel;
                break;
            }
            default:
            {
                NEVER_REACHED;
            }
        }
// LCOV_EXCL_START
        return NULL;  //This is included to appease compilers which analysis the NEVER_REACHED branch as not returning a pointer
// LCOV_EXCL_STOP
    }
};

#endif // LABELBASEDCONTRACTIONCELLFACTORY_HPP_
