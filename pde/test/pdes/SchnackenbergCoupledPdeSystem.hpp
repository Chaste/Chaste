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

#ifndef SCHNACKENBERGCOUPLEDPDESYSTEM_HPP_
#define SCHNACKENBERGCOUPLEDPDESYSTEM_HPP_

#include "AbstractLinearParabolicPdeSystemForCoupledOdeSystem.hpp"
#include "ChastePoint.hpp"
#include "Element.hpp"

/**
 * Two coupled PDEs defining the Schnackenberg reaction-diffusion system
 *
 * u_t = D1*del^2 u + kappa1 - kappa_1*u + kappa3*u^2*v,
 * v_t = D2*del^2 u + kappa2 - kappa3*u^2*v.
 */
template<unsigned SPACE_DIM>
class SchnackenbergCoupledPdeSystem : public AbstractLinearParabolicPdeSystemForCoupledOdeSystem<SPACE_DIM, SPACE_DIM, 2>
{
private:

    double mD1;      /**< Parameter D1 in the Schnackenberg system. */
    double mD2;      /**< Parameter D2 in the Schnackenberg system. */
    double mKappa1;  /**< Parameter kappa1 in the Schnackenberg system. */
    double mKappa_1; /**< Parameter kappa_1 in the Schnackenberg system. */
    double mKappa2;  /**< Parameter kappa2 in the Schnackenberg system. */
    double mKappa3;  /**< Parameter kappa3 in the Schnackenberg system. */

public:

    SchnackenbergCoupledPdeSystem(double d1=1.0,
                                  double d2=1.0,
                                  double kappa1=1.0,
                                  double kappa_1=1.0,
                                  double kappa2=1.0,
                                  double kappa3=1.0)
        : AbstractLinearParabolicPdeSystemForCoupledOdeSystem<SPACE_DIM, SPACE_DIM, 2>(),
          mD1(d1),
          mD2(d2),
          mKappa1(kappa1),
          mKappa_1(kappa_1),
          mKappa2(kappa2),
          mKappa3(kappa3)
    {
    }

    double ComputeDuDtCoefficientFunction(const ChastePoint<SPACE_DIM>& rX, unsigned index)
    {
        return 1.0;
    }

    double ComputeSourceTerm(const ChastePoint<SPACE_DIM>& rX, c_vector<double,2>& rU, std::vector<double>& rOdeSolution, unsigned pdeIndex)
    {
        assert(pdeIndex == 0 || pdeIndex == 1);

        double source_term;
        if (pdeIndex == 0)
        {
            source_term = mKappa1 - mKappa_1*rU(0) + mKappa3*rU(1)*rU(0)*rU(0);
        }
        else // pdeIndex == 1
        {
            source_term = mKappa2 - mKappa3*rU(1)*rU(0)*rU(0);
        }
        return source_term;
    }

    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& rX, unsigned pdeIndex, Element<SPACE_DIM,SPACE_DIM>* pElement=NULL)
    {
        assert(pdeIndex == 0 || pdeIndex == 1);

        c_matrix<double, SPACE_DIM, SPACE_DIM> diffusion_term;
        if (pdeIndex == 0)
        {
            diffusion_term = mD1*identity_matrix<double>(SPACE_DIM);
        }
        else // pdeIndex == 1
        {
            diffusion_term = mD2*identity_matrix<double>(SPACE_DIM);
        }
        return diffusion_term;
    }
};

#endif /*SCHNACKENBERGCOUPLEDPDESYSTEM_HPP_*/
