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
     * @return the first derivative dW/dI1.
     *
     * \todo The name of this method should not include underscores.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    double Get_dW_dI1(double I1, double I2);

    /**
     * @return the first derivative dW/dI2.
     *
     * \todo The name of this method should not include underscores.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    double Get_dW_dI2(double I1, double I2);

    /**
     * @return the second derivative d^2W/dI1^2.
     *
     * \todo The name of this method should not include underscores.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    double Get_d2W_dI1(double I1, double I2);

    /**
     * @return the second derivative d^2W/dI2^2.
     *
     * \todo The name of this method should not include underscores.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    double Get_d2W_dI2(double I1, double I2);

    /**
     * @return the second derivative d^2W/dI1dI2.
     *
     * \todo The name of this method should not include underscores.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    double Get_d2W_dI1I2(double I1, double I2);

    /**
     * @return the parameter alpha_{ij}.
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
     * @return the matrix alpha
     */
    static std::vector<std::vector<double> > GetZeroedAlpha(unsigned n);
};

#endif /*POLYNOMIALMATERIALLAW3D_HPP_*/
