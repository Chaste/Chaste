/*

Copyright (C) University of Oxford, 2005-2012

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

#ifndef _ABSTRACTLINEARELLIPTICPDE_HPP_
#define _ABSTRACTLINEARELLIPTICPDE_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"

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
class AbstractLinearEllipticPde
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
    }

public:

    /**
     * Constructor.
     */
    AbstractLinearEllipticPde()
    {}

    /**
     * Destructor.
     */
    virtual ~AbstractLinearEllipticPde()
    {}

    /**
     * Compute the constant in u part of the source term, i.e g(x) in
     * Div(D Grad u)  +  f(x)u + g(x) = 0, at a given point.
     *
     * @param rX The point in space
     * @param pElement The element
     */
    virtual double ComputeConstantInUSourceTerm(const ChastePoint<SPACE_DIM>& rX, Element<ELEMENT_DIM,SPACE_DIM>* pElement)=0;

    /**
     * Compute the coefficient of u in the linear part of the source term, i.e f(x) in
     * Div(D Grad u)  +  f(x)u + g(x) = 0, at a given point in space.
     *
     * @param rX The point in space
     * @param pElement
     */
    virtual double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<SPACE_DIM>& rX,
                                                     Element<ELEMENT_DIM,SPACE_DIM>* pElement)=0;

    /**
     * Compute the diffusion term at a given point. The diffusion tensor should be symmetric and positive definite
     *
     * @param rX The point in space at which the diffusion term is computed.
     * @return A matrix.
     */
    virtual c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& rX)=0;

    /**
     * Compute the constant in u part of the source term, i.e g(x) in
     * Div(D Grad u)  +  f(x)u + g(x) = 0, at a given node.
     *
     * @param rNode the node
     */
    virtual double ComputeConstantInUSourceTermAtNode(const Node<SPACE_DIM>& rNode);

    /**
     * Compute the coefficient of u in the linear part of the source term, i.e f(x) in
     * Div(D Grad u)  +  f(x)u + g(x) = 0, at a given node.
     *
     * @param rNode the node
     */
    virtual double ComputeLinearInUCoeffInSourceTermAtNode(const Node<SPACE_DIM>& rNode);
};

TEMPLATED_CLASS_IS_ABSTRACT_2_UNSIGNED(AbstractLinearEllipticPde)

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::ComputeConstantInUSourceTermAtNode(const Node<SPACE_DIM>& rNode)
{
    return ComputeConstantInUSourceTerm(rNode.GetPoint(), NULL);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>::ComputeLinearInUCoeffInSourceTermAtNode(const Node<SPACE_DIM>& rNode)
{
    return ComputeLinearInUCoeffInSourceTerm(rNode.GetPoint(), NULL);
}

#endif //_ABSTRACTLINEARELLIPTICPDE_HPP_
