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

#ifndef _VARYINGDIFFUSIONANDSOURCETERMPDE_
#define _VARYINGDIFFUSIONANDSOURCETERMPDE_

#include "AbstractLinearEllipticPde.hpp"
#include "ChastePoint.hpp"
#include <cmath>

/**
 * A more complex linear elliptic PDE used in tests. The source and diffusion terms
 * depend on x.
 */
template <int SPACE_DIM>
class VaryingDiffusionAndSourceTermPde : public AbstractLinearEllipticPde<SPACE_DIM,SPACE_DIM>
{
private:
    double DistanceFromOrigin(const ChastePoint<SPACE_DIM>& x)
    {
        double sum=0;
        for (int i=0; i<SPACE_DIM; i++)
        {
            sum += x[i]*x[i];
        }
        return sqrt(sum);
    }

public:
    double ComputeConstantInUSourceTerm(const ChastePoint<SPACE_DIM>& rX, Element<SPACE_DIM,SPACE_DIM>*)
    {
        return pow(DistanceFromOrigin(rX),3);
    }

    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<SPACE_DIM>& , Element<SPACE_DIM,SPACE_DIM>*)
    {
        return 0.0;
    }

    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& x)
    {
        return pow(DistanceFromOrigin(x),2)*identity_matrix<double>(SPACE_DIM);
    }
};

#endif //_VARYINGDIFFUSIONANDSOURCETERMPDE_
