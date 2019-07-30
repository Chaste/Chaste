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

#ifndef _ABSTRACTLINEARELLIPTICPDE_HPP_
#define _ABSTRACTLINEARELLIPTICPDE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractLinearPde.hpp"
#include "UblasCustomFunctions.hpp"
#include "ChastePoint.hpp"
#include "Node.hpp"
#include "Element.hpp"
#include <petscvec.h>

/**
 * AbstractLinearEllipticPde class.
 *
 * A general PDE of the form:
 * 0 =   Grad.(DiffusionTerm(x)*Grad(u))
 *     + ComputeConstantInUSourceTerm(x)
 *     + ComputeLinearInUCoeffInSourceTerm(x, u)
 *
 * Parabolic PDEs are be derived from this (AbstractLinearParabolicPde)
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractLinearEllipticPde : public AbstractLinearPde<ELEMENT_DIM, SPACE_DIM>
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
    AbstractLinearEllipticPde()
        : AbstractLinearPde<ELEMENT_DIM, SPACE_DIM>()
    {}

    /**
     * Destructor.
     */
    virtual ~AbstractLinearEllipticPde()
    {}

    /**
     * @return computed constant in u part of the source term, i.e g(x) in
     * Div(D Grad u)  +  f(x)u + g(x) = 0, at a given point.
     *
     * @param rX The point in space
     * @param pElement The element
     */
    virtual double ComputeConstantInUSourceTerm(const ChastePoint<SPACE_DIM>& rX,
                                                Element<ELEMENT_DIM,SPACE_DIM>* pElement)=0;

    /**
     * @return computed coefficient of u in the linear part of the source term, i.e f(x) in
     * Div(D Grad u)  +  f(x)u + g(x) = 0, at a given point in space.
     *
     * @param rX The point in space
     * @param pElement
     */
    virtual double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<SPACE_DIM>& rX,
                                                     Element<ELEMENT_DIM,SPACE_DIM>* pElement)=0;

    /**
     * @return computed diffusion term at a given point. The diffusion tensor should be symmetric and positive definite
     *
     * @param rX The point in space at which the diffusion term is computed.
     * @return A matrix.
     */
    virtual c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& rX)=0;

    /**
     * @return computed constant in u part of the source term, i.e g(x) in
     * Div(D Grad u)  +  f(x)u + g(x) = 0, at a given node.
     *
     * @param rNode the node
     */
    virtual double ComputeConstantInUSourceTermAtNode(const Node<SPACE_DIM>& rNode);

    /**
     * @return computed coefficient of u in the linear part of the source term, i.e f(x) in
     * Div(D Grad u)  +  f(x)u + g(x) = 0, at a given node.
     *
     * @param rNode the node
     */
    virtual double ComputeLinearInUCoeffInSourceTermAtNode(const Node<SPACE_DIM>& rNode);
};

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::ComputeConstantInUSourceTermAtNode(const Node<SPACE_DIM>& rNode)
{
    return ComputeConstantInUSourceTerm(rNode.GetPoint(), nullptr);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::ComputeLinearInUCoeffInSourceTermAtNode(const Node<SPACE_DIM>& rNode)
{
    return ComputeLinearInUCoeffInSourceTerm(rNode.GetPoint(), nullptr);
}

#endif //_ABSTRACTLINEARELLIPTICPDE_HPP_
