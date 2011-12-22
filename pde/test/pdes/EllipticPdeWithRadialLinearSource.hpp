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

#ifndef ELLIPTICPDEWITHRADIALLINEARSOURCE_HPP_
#define ELLIPTICPDEWITHRADIALLINEARSOURCE_HPP_

/**
 * The pde u_xx+u_yy - (x^2+y^2)u = 0, just for one particular test. This has the solution
 * exp(xy)
 */
class EllipticPdeWithRadialLinearSource :public AbstractLinearEllipticPde<2,2>
{
public:

    double ComputeConstantInUSourceTerm(const ChastePoint<2>&, Element<2,2>*)
    {
        return 0.0;
    }

    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<2>& rX, Element<2,2>*)
    {
        return -(rX[0]*rX[0] + rX[1]*rX[1]);
    }

    c_matrix<double, 2, 2> ComputeDiffusionTerm(const ChastePoint<2>& )
    {
        return identity_matrix<double>(2);
    }
};

#endif /*ELLIPTICPDEWITHRADIALLINEARSOURCE_HPP_*/
