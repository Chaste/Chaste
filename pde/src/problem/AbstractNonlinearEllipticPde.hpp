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
#ifndef _ABSTRACTNONLINEARELLIPTICPDE_HPP_
#define _ABSTRACTNONLINEARELLIPTICPDE_HPP_

#include "UblasCustomFunctions.hpp"
#include "ChastePoint.hpp"


/**
 * AbstractNonlinearEllipticPde class.
 *
 * A simple elliptic PDE in 1 unknown with nonlinear diffusion term as
 * well as nonlinear source term:
 *
 *  0 = Grad.(DiffusionTerm(x,u)*Grad(u))+LinearSourceTerm(x)+NonlinearSourceTerm(x, u)
 *
 */
template <unsigned SPACE_DIM>
class AbstractNonlinearEllipticPde
{
public:

    /**
     * Compute linear source term.
     *
     * @param rX the point in space at which the linear source term is computed
     */
    virtual double ComputeLinearSourceTerm(const ChastePoint<SPACE_DIM>& rX)=0;

    /**
     * Compute nonlinear source term.
     *
     * @param rX the point in space at which the nonlinear source term is computed
     * @param u the value of the dependent variable at the point
     */
    virtual double ComputeNonlinearSourceTerm(const ChastePoint<SPACE_DIM>& rX, double u)=0;

    /**
     * Compute diffusion term. The diffusion tensor should be symmetric and positive definite.
     *
     * @param rX the point in space at which the diffusion term is computed.
     * @param u the value of the dependent variable at the point
     */
    virtual c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& rX, double u)=0;

    /**
     * Compute derivative of diffusion term.
     *
     * @param rX the point in space at which the diffusion term is computed.
     * @param u the value of the dependent variable at the point
     */
    virtual c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTermPrime(const ChastePoint<SPACE_DIM>& rX, double u)=0;

    /**
     * Compute derivative of nonlinear source term.
     *
     * @param rX the point in space at which the nonlinear source term is computed
     * @param u the value of the dependent variable at the point
     */
    virtual double ComputeNonlinearSourceTermPrime(const ChastePoint<SPACE_DIM>& rX, double u)=0;

    /**
     * Destructor.
     */
    virtual ~AbstractNonlinearEllipticPde()
    {}
};

#endif //_ABSTRACTNONLINEARELLIPTICPDE_HPP_
