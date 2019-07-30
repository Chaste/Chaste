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
#ifndef _ABSTRACTLINEARPARABOLICPDE_HPP_
#define _ABSTRACTLINEARPARABOLICPDE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractLinearPde.hpp"
#include "UblasCustomFunctions.hpp"
#include "ChastePoint.hpp"
#include "Node.hpp"
#include "Element.hpp"
#include <petscvec.h>

/**
 * AbstractLinearParabolicPde class.
 *
 * A general PDE of the form:
 * c(x) du/dt = Grad.(DiffusionTerm(x)*Grad(u))+LinearSourceTerm(x)+NonlinearSourceTerm(x, u)
 *
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class AbstractLinearParabolicPde : public AbstractLinearPde<ELEMENT_DIM, SPACE_DIM>
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the PDE object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractLinearPde<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    /**
     * Constructor.
     */
    AbstractLinearParabolicPde()
        : AbstractLinearPde<ELEMENT_DIM, SPACE_DIM>()
    {}

    /**
     * Destructor.
     */
    virtual ~AbstractLinearParabolicPde()
    {}

    /**
     * @return the function c(x) in "c(x) du/dt = Grad.(DiffusionTerm(x)*Grad(u))+LinearSourceTerm(x)+NonlinearSourceTerm(x, u)"
     *
     * @param rX the point in space at which the function c is computed
     */
    virtual double ComputeDuDtCoefficientFunction(const ChastePoint<SPACE_DIM>& rX)=0;

    /**
     * @return computed source term.
     *
     * @param rX the point in space at which the nonlinear source term is computed
     * @param u the value of the dependent variable at the point
     * @param pElement the element that we are inside
     */
    virtual double ComputeSourceTerm(const ChastePoint<SPACE_DIM>& rX,
                                     double u,
                                     Element<ELEMENT_DIM,SPACE_DIM>* pElement=nullptr)=0;

    /**
     * @return computed source term at a node.
     *
     * @param rNode the node at which the nonlinear source term is computed
     * @param u the value of the dependent variable at the node
     */
    virtual double ComputeSourceTermAtNode(const Node<SPACE_DIM>& rNode,
                                           double u);


    /**
     * @return computed diffusion term. The diffusion tensor should be symmetric and positive definite.
     *
     * @param rX The point in space at which the diffusion term is computed.
     * @param pElement The mesh element that x is contained in (optional).
     * @return A matrix.
     */
    virtual c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& rX,
                                                                        Element<ELEMENT_DIM,SPACE_DIM>* pElement=NULL)=0;
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractLinearParabolicPde<ELEMENT_DIM, SPACE_DIM>::ComputeSourceTermAtNode(const Node<SPACE_DIM>& rNode,
                                                                                   double u)
{
    return ComputeSourceTerm(rNode.GetPoint(), u);
}

#endif //_ABSTRACTLINEARPARABOLICPDE_HPP_
