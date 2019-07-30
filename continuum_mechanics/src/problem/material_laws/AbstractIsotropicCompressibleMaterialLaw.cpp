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

#include "AbstractIsotropicCompressibleMaterialLaw.hpp"

template<unsigned DIM>
AbstractIsotropicCompressibleMaterialLaw<DIM>::~AbstractIsotropicCompressibleMaterialLaw()
{
}

template<unsigned DIM>
void AbstractIsotropicCompressibleMaterialLaw<DIM>::ComputeStressAndStressDerivative(c_matrix<double,DIM,DIM>& rC,
                                                                                     c_matrix<double,DIM,DIM>& rInvC,
                                                                                     double                    pressure,
                                                                                     c_matrix<double,DIM,DIM>& rT,
                                                                                     FourthOrderTensor<DIM,DIM,DIM,DIM>&   rDTdE,
                                                                                     bool                      computeDTdE)
{
    /*
     * This is covered, but gcov doesn't see this as being covered
     * because it runs at compile time checking templates I think.
     */
    assert((DIM==2) || (DIM==3)); // LCOV_EXCL_LINE
    assert(pressure==0.0);

    static c_matrix<double,DIM,DIM> identity = identity_matrix<double>(DIM);

    double I1 = Trace(rC);
    double I2 = SecondInvariant(rC);
    double I3 = Determinant(rC);

    static c_matrix<double,DIM,DIM> dI2dC;
    dI2dC = I1*identity - rC;              // MUST be on separate line to above!

    double w1 = Get_dW_dI1(I1,I2,I3);
    double w2 = Get_dW_dI2(I1,I2,I3);
    double w3 = Get_dW_dI3(I1,I2,I3);


    // Compute stress:  **** See FiniteElementImplementations document. ****
    //
    //  T = dW_dE
    //    = 2 dW_dC
    //    = 2 (  w1 dI1/dC   +  w2 dI2/dC      +   w3 dI3/dC )
    //    = 2 (  w1 I        +  w2 (I1*I - C)  +   w3 I3 inv(C) )
    //
    //  where w1 = dW/dI1, etc
    //
    rT = 2*w1*identity + 2*w3*I3*rInvC;
    if (DIM==3)
    {
        rT += 2*w2*dI2dC;
    }

    // Compute stress derivative if required:  **** See FiniteElementImplementations document. ****
    //
    // The stress derivative dT_{MN}/dE_{PQ} is
    //
    //
    //  dT_dE = 2 dT_dC
    //        = 4  d/dC ( w1 I  +  w2 (I1*I - C)  +   w3 I3 inv(C) )
    //  so (in the following ** represents outer product):
    //  (1/4) dT_dE =        w11 I**I          +    w12 I**(I1*I-C)           +     w13 I**inv(C)
    //                  +    w21 (I1*I-C)**I   +    w22 (I1*I-C)**(I1*I-C)    +     w23 (I1*I-C)**inv(C)           +   w2 (I**I - dC/dC)
    //                  +    w31 I3 inv(C)**I  +    w32 I3 inv(C)**(I1*I-C)   +  (w33 I3 + w3) inv(C)**inv(C)      +   w3 d(invC)/dC
    //
    //  Here, I**I represents the tensor A[M][N][P][Q] = (M==N)*(P==Q) // ie delta(M,N)delta(P,Q),   etc
    //

    if (computeDTdE)
    {
        double  w11    = Get_d2W_dI1(I1,I2,I3);
        double  w22    = Get_d2W_dI2(I1,I2,I3);
        double  w33    = Get_d2W_dI3(I1,I2,I3);

        double  w23  = Get_d2W_dI2I3(I1,I2,I3);
        double  w13  = Get_d2W_dI1I3(I1,I2,I3);
        double  w12  = Get_d2W_dI1I2(I1,I2,I3);

        for (unsigned M=0; M<DIM; M++)
        {
            for (unsigned N=0; N<DIM; N++)
            {
                for (unsigned P=0; P<DIM; P++)
                {
                    for (unsigned Q=0; Q<DIM; Q++)
                    {
                        rDTdE(M,N,P,Q) =   4 * w11  * (M==N) * (P==Q)
                                         + 4 * w13  * I3 * ( (M==N) * rInvC(P,Q)  +  rInvC(M,N)*(P==Q) )  // the w13 and w31 terms
                                         + 4 * (w33*I3 + w3) * I3 * rInvC(M,N) * rInvC(P,Q)
                                         - 4 * w3 * I3 * rInvC(M,P) * rInvC(Q,N);

                        if (DIM==3)
                        {
                            rDTdE(M,N,P,Q) +=   4 * w22  * dI2dC(M,N) * dI2dC(P,Q)
                                              + 4 * w12  * ((M==N)*dI2dC(P,Q) + (P==Q)*dI2dC(M,N))          // the w12 and w21 terms
                                              + 4 * w23 * I3 * ( dI2dC(M,N)*rInvC(P,Q) + rInvC(M,N)*dI2dC(P,Q)) // the w23 and w32 terms
                                              + 4 * w2   * ((M==N)*(P==Q) - (M==P)*(N==Q));
                        }
                    }
                }
            }
        }
    }
}

// Explicit instantiation
template class AbstractIsotropicCompressibleMaterialLaw<2>;
template class AbstractIsotropicCompressibleMaterialLaw<3>;
