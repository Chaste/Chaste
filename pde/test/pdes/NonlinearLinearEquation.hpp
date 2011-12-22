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

#ifndef _NONLINEARLINEAREQUATION_HPP_
#define _NONLINEARLINEAREQUATION_HPP_

#include "AbstractNonlinearEllipticPde.hpp"

/**
 * Linear poisson equation, using the nonlinear class format.
 */
template <int SPACE_DIM>
class NonlinearLinearEquation : public AbstractNonlinearEllipticPde<SPACE_DIM>
{
public:

    double ComputeLinearSourceTerm(const ChastePoint<SPACE_DIM>& )
    {
        return 1.0;
    }

    double ComputeNonlinearSourceTerm(const ChastePoint<SPACE_DIM>& , double )
    {
        return 0.0;
    }

    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& , double u)
    {
        return identity_matrix<double>(SPACE_DIM);
    }

    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTermPrime(const ChastePoint<SPACE_DIM>& , double )
    {
        return identity_matrix<double>(SPACE_DIM)*0.0;
    }

    double ComputeNonlinearSourceTermPrime(const ChastePoint<SPACE_DIM>& , double )
    {
        return 0.0;
    }
};

#endif //_NONLINEARLINEAREQUATION_HPP_
