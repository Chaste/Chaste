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
     * Compute the function c_i(x).
     *
     * @param rX the point x at which the function c_i is computed
     * @param pdeIndex the index of the PDE, denoted by i above
     */
    virtual double ComputeDuDtCoefficientFunction(const ChastePoint<SPACE_DIM>& rX, unsigned pdeIndex)=0;

    /**
     * Compute the source term f_i(x, u_1, u_2, ..., u_p) at a point in space.
     *
     * @param rX the point x at which the nonlinear source term is computed
     * @param rU the vector of dependent variables (u_1, u_2, ..., u_p) at the point x
     * @param rOdeSolution the ODE system state vector (v_1, ..., v_q) at the point x (if an ODE system is present)
     * @param pdeIndex the index of the PDE, denoted by i above
     */
    virtual double ComputeSourceTerm(const ChastePoint<SPACE_DIM>& rX, c_vector<double,PROBLEM_DIM>& rU, std::vector<double>& rOdeSolution, unsigned pdeIndex)=0;

    /**
     * Compute source term f_i(x, u_1, u_2, ..., u_p) at a node.
     *
     * @param rNode the node at which the nonlinear source term f_i is computed
     * @param rU the vector of dependent variables (u_1, u_2, ..., u_p) at the node
     * @param rOdeSolution the ODE system state vector (v_1, ..., v_q) at the node (if an ODE system is present)
     * @param pdeIndex the index of the PDE, denoted by i above
     */
    virtual double ComputeSourceTermAtNode(const Node<SPACE_DIM>& rNode, c_vector<double,PROBLEM_DIM>& rU, std::vector<double>& rOdeSolution, unsigned pdeIndex);

    /**
     * Compute diffusion term D_i(x) at a point in space. The diffusion tensor should be symmetric and positive definite.
     *
     * @param rX The point x at which the diffusion term D_i is computed
     * @param pdeIndex the index of the PDE, denoted by i above
     * @param pElement The mesh element that x is contained in (optional)
     */
    virtual c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& rX, unsigned pdeIndex, Element<ELEMENT_DIM,SPACE_DIM>* pElement=NULL)=0;

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
