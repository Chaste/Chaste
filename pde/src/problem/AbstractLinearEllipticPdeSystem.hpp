/*

Copyright (c) 2005-2017, University of Oxford.
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

#ifndef _ABSTRACTLINEARELLIPTICPDESYSTEM_HPP_
#define _ABSTRACTLINEARELLIPTICPDESYSTEM_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractLinearPdeSystem.hpp"
#include "UblasCustomFunctions.hpp"
#include "ChastePoint.hpp"
#include "Node.hpp"
#include "Element.hpp"
#include <petscvec.h>

/**
 * A system of (uncoupled) elliptic PDEs with linear source terms:
 *
 * 0 = div (D_i(x) grad (u_i)) + f_i(x)u_i + g_i(x),  i=0,...,p-1.
 *
 * Here p denotes the size of the PDE system.
 *
 * Such systems may be solved using LinearEllipticPdeSystemSolver.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM=1>
class AbstractLinearEllipticPdeSystem : public AbstractLinearPdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
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
        archive & boost::serialization::base_object<AbstractLinearPdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM> >(*this);
    }

public:

    /**
     * Constructor.
     */
    AbstractLinearEllipticPdeSystem()
        : AbstractLinearPdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>()
    {}

    /**
     * Destructor.
     */
    virtual ~AbstractLinearEllipticPdeSystem()
    {}

    /**
     * @return computed contribution to source term that is constant in u, g_i(x), at a point in space.
     *
     * @param rX the point x at which the source term contribution is computed
     * @param pdeIndex the index of the PDE, denoted by i above
     * @param pElement The mesh element that x is contained in
     */
    virtual double ComputeConstantInUSourceTerm(const ChastePoint<SPACE_DIM>& rX,
                                                unsigned pdeIndex,
                                                Element<ELEMENT_DIM,SPACE_DIM>* pElement)=0;

    /**
     * @return computed contribution to source term that is linear in u, f_i(x), at a point in space.
     *
     * @param rX The point in space
     * @param pdeIndex the index of the PDE, denoted by i above
     * @param pElement The mesh element that x is contained in
     */
    virtual double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<SPACE_DIM>& rX,
                                                     unsigned pdeIndex,
                                                     Element<ELEMENT_DIM,SPACE_DIM>* pElement)=0;

    /**
     * @return computed diffusion term D_i(x) at a point in space. The diffusion tensor should be symmetric and positive definite.
     *
     * @param rX The point x at which the diffusion term D_i is computed
     * @param pdeIndex the index of the PDE, denoted by i above
     */
    virtual c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& rX,
                                                                        unsigned pdeIndex)=0;

    /**
     * @return computed contribution to source term that is constant in u, g_i(x), at a node.
     *
     * @param rNode the node at which the source term contribution is computed
     * @param pdeIndex the index of the PDE, denoted by i above
     */
    virtual double ComputeConstantInUSourceTermAtNode(const Node<SPACE_DIM>& rNode,
                                                      unsigned pdeIndex);

    /**
     * @return computed contribution to source term that is linear in u, f_i(x), at a node.
     *
     * @param rNode the node at which the source term contribution is computed
     * @param pdeIndex the index of the PDE, denoted by i above
     */
    virtual double ComputeLinearInUCoeffInSourceTermAtNode(const Node<SPACE_DIM>& rNode,
                                                           unsigned pdeIndex);
};

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
double AbstractLinearEllipticPdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ComputeConstantInUSourceTermAtNode(const Node<SPACE_DIM>& rNode, unsigned pdeIndex)
{
    return ComputeConstantInUSourceTerm(rNode.GetPoint(), pdeIndex, nullptr);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
double AbstractLinearEllipticPdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ComputeLinearInUCoeffInSourceTermAtNode(const Node<SPACE_DIM>& rNode, unsigned pdeIndex)
{
    return ComputeLinearInUCoeffInSourceTerm(rNode.GetPoint(), pdeIndex, nullptr);
}

#endif //_ABSTRACTLINEARELLIPTICPDESYSTEM_HPP_
