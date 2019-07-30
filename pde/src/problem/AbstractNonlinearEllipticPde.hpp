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
     * @return computed linear source term.
     *
     * @param rX the point in space at which the linear source term is computed
     */
    virtual double ComputeLinearSourceTerm(const ChastePoint<SPACE_DIM>& rX)=0;

    /**
     * @return computed nonlinear source term.
     *
     * @param rX the point in space at which the nonlinear source term is computed
     * @param u the value of the dependent variable at the point
     */
    virtual double ComputeNonlinearSourceTerm(const ChastePoint<SPACE_DIM>& rX, double u)=0;

    /**
     * @return computed diffusion term. The diffusion tensor should be symmetric and positive definite.
     *
     * @param rX the point in space at which the diffusion term is computed.
     * @param u the value of the dependent variable at the point
     */
    virtual c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& rX, double u)=0;

    /**
     * @return computed derivative of diffusion term.
     *
     * @param rX the point in space at which the diffusion term is computed.
     * @param u the value of the dependent variable at the point
     */
    virtual c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTermPrime(const ChastePoint<SPACE_DIM>& rX, double u)=0;

    /**
     * @return computed derivative of nonlinear source term.
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
