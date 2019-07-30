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

#ifndef ABSTRACTLINEARPARABOLICPDESYSTEMFORCOUPLEDODESYSTEM_HPP_
#define ABSTRACTLINEARPARABOLICPDESYSTEMFORCOUPLEDODESYSTEM_HPP_

#include "UblasCustomFunctions.hpp"
#include "ChastePoint.hpp"
#include "Node.hpp"
#include "Element.hpp"
#include <petscvec.h>

/**
 * AbstractLinearParabolicPdeSystemForCoupledOdeSystem class.
 *
 * A system of parabolic PDEs, which may be coupled via their source terms:
 *
 * d/dt (u_i) = div (D(x) grad (u_i)) + f_i (x, u_0, ..., u_{p-1}, v_0, ..., v_{q-1}),  i=0,...,p-1.
 *
 * Here p denotes the size of the PDE system and each source term f_i may be nonlinear.
 * The variables v_0, ..., v_{q-1} are assumed to satisfy a coupled ODE system of the form
 *
 * d/dt (v_j) = g_j(x, u_0, ..., u_{p-1}, v_0, ..., v_{q-1}),  j=0,...,q-1.
 *
 * Such systems may be solved using LinearParabolicPdeSystemWithCoupledOdeSystemSolver.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM, unsigned PROBLEM_DIM=1>
class AbstractLinearParabolicPdeSystemForCoupledOdeSystem
{
public:
    /**
     * @return computed function c_i(x).
     *
     * @param rX the point x at which the function c_i is computed
     * @param pdeIndex the index of the PDE, denoted by i above
     */
    virtual double ComputeDuDtCoefficientFunction(const ChastePoint<SPACE_DIM>& rX, unsigned pdeIndex)=0;

    /**
     * @return computed source term f_i(x, u_1, u_2, ..., u_p) at a point in space.
     *
     * @param rX the point x at which the nonlinear source term is computed
     * @param rU the vector of dependent variables (u_1, u_2, ..., u_p) at the point x
     * @param rOdeSolution the ODE system state vector (v_1, ..., v_q) at the point x (if an ODE system is present)
     * @param pdeIndex the index of the PDE, denoted by i above
     */
    virtual double ComputeSourceTerm(const ChastePoint<SPACE_DIM>& rX,
                                     c_vector<double,PROBLEM_DIM>& rU,
                                     std::vector<double>& rOdeSolution,
                                     unsigned pdeIndex)=0;

    /**
     * @return computed source term f_i(x, u_1, u_2, ..., u_p) at a node.
     *
     * @param rNode the node at which the nonlinear source term f_i is computed
     * @param rU the vector of dependent variables (u_1, u_2, ..., u_p) at the node
     * @param rOdeSolution the ODE system state vector (v_1, ..., v_q) at the node (if an ODE system is present)
     * @param pdeIndex the index of the PDE, denoted by i above
     */
    virtual double ComputeSourceTermAtNode(const Node<SPACE_DIM>& rNode,
                                           c_vector<double,PROBLEM_DIM>& rU,
                                           std::vector<double>& rOdeSolution,
                                           unsigned pdeIndex);

    /**
     * @return computed diffusion term D_i(x) at a point in space. The diffusion tensor should be symmetric and positive definite.
     *
     * @param rX The point x at which the diffusion term D_i is computed
     * @param pdeIndex the index of the PDE, denoted by i above
     * @param pElement The mesh element that x is contained in (optional)
     */
    virtual c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& rX,
                                                                        unsigned pdeIndex,
                                                                        Element<ELEMENT_DIM,SPACE_DIM>* pElement=NULL)=0;

    /**
     * Destructor.
     */
    virtual ~AbstractLinearParabolicPdeSystemForCoupledOdeSystem()
    {}
};

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
double AbstractLinearParabolicPdeSystemForCoupledOdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::ComputeSourceTermAtNode(const Node<SPACE_DIM>& rNode, c_vector<double,PROBLEM_DIM>& rU, std::vector<double>& rOdeSolution, unsigned pdeIndex)
{
    return ComputeSourceTerm(rNode.GetPoint(), rU, rOdeSolution, pdeIndex);
}

#endif /*ABSTRACTLINEARPARABOLICPDESYSTEMFORCOUPLEDODESYSTEM_HPP_*/
