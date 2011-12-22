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


#ifndef POLYNOMIALMATERIALLAW3D_HPP_
#define POLYNOMIALMATERIALLAW3D_HPP_

#include "AbstractIsotropicIncompressibleMaterialLaw.hpp"


/**
 *  PolynomialMaterialLaw3d
 *
 *  An incompressible, isotropic, hyperelastic material law with a polynomial form
 *
 *  W(I_1,I_2)  =  Sigma_{0<p+q<=N}  alpha_{pq} (I_1-3)^p (I_2-3)^q   -  (pressure/2) C^{-1}
 *
 *  For example, if N=1, this reduces to the Mooney Rivlin law
 *     W(I_1,I_2)  =  alpha_{10} (I_1-3) +  alpha_{01} (I_2-3)   -  (pressure/2) C^{-1}
 *  ie the matrix alpha has the form
 *  [ 0  c1 ]
 *  [ c2  0 ]
 *  where c1 and c2 is the usual notation for the Mooney-Rivlin constants
 *
 *  The polynomial is specified by passing in N and the matrix (actually a std::vector
 *  of std::vector<double>s) alpha. alpha should be of size N+1 by N+1, with the bottom
 *  right hand block (ie the components such that p+q>N) all zero. alpha[0][0] should
 *  really also be 0, but, being since alpha[0][0] (I1_3)^0 (I2-3)^0 is a constant
 *  and disappears when the strain energy W is differentiated to obtain the stress, it is
 *  not used. An exception is thrown if alpha[p][q]!=0 for p+q > N though.
 */
class PolynomialMaterialLaw3d : public AbstractIsotropicIncompressibleMaterialLaw<3>
{
private:

    /** Parameter N. */
    unsigned mN;

    /** Matrix of parameters alpha. */
    std::vector< std::vector<double> > mAlpha;

public:

    /**
     * Get the first derivative dW/dI1.
     *
     * \todo The name of this method should not include underscores.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    double Get_dW_dI1(double I1, double I2);

    /**
     * Get the first derivative dW/dI2.
     *
     * \todo The name of this method should not include underscores.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    double Get_dW_dI2(double I1, double I2);

    /**
     * Get the second derivative d^2W/dI1^2.
     *
     * \todo The name of this method should not include underscores.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    double Get_d2W_dI1(double I1, double I2);

    /**
     * Get the second derivative d^2W/dI2^2.
     *
     * \todo The name of this method should not include underscores.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    double Get_d2W_dI2(double I1, double I2);

    /**
     * Get the second derivative d^2W/dI1dI2.
     *
     * \todo The name of this method should not include underscores.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    double Get_d2W_dI1I2(double I1, double I2);

    /**
     * Get the parameter alpha_{ij}.
     *
     * @param i index i
     * @param j index j
     */
    double GetAlpha(unsigned i, unsigned j);

public:

    /**
     * Constructor.
     *
     * @param n the parameter n
     * @param alpha the matrix of parameters alpha
     */
    PolynomialMaterialLaw3d(unsigned n, std::vector<std::vector<double> > alpha);

    /**
     * Resize the matrix alpha to be of size (n+1)*(n+1) and zero all entries.
     *
     * @param n the parameter n
     */
    static std::vector<std::vector<double> > GetZeroedAlpha(unsigned n);
};

#endif /*POLYNOMIALMATERIALLAW3D_HPP_*/
