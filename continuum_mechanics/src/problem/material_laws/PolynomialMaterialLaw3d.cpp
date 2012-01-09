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

#include "PolynomialMaterialLaw3d.hpp"

double PolynomialMaterialLaw3d::Get_dW_dI1(double I1, double I2)
{
    double ret = 0.0;
    // notes: use ints not unsigned as doing p-1
    // (except indexing from p=1 because multiplying by p, but
    // still safer to use ints)
    for (int p=1; p<=(int)mN; p++)
    {
        for (int q=0; q<=(int)mN-p; q++)
        {
            ret += mAlpha[p][q] * p * pow(I1-3,p-1) * pow(I2-3,q);
        }
    }

    return ret;
}

double PolynomialMaterialLaw3d::Get_dW_dI2(double I1, double I2)
{
    double ret = 0.0;
    // notes: use ints not unsigned as doing q-1
    // (except indexing from q=1 because multiplying by q, but
    // still safer to use ints)
    for (int p=0; p<=(int)mN; p++)
    {
        for (int q=1; q<=(int)mN-p; q++)
        {
            ret += mAlpha[p][q] * q * pow(I1-3,p) * pow(I2-3,q-1);
        }
    }
    return ret;
}

double PolynomialMaterialLaw3d::Get_d2W_dI1(double I1, double I2)
{
    double ret = 0.0;

    // notes: use ints not unsigned as doing p-1
    // (except indexing from p=2 because multiplying by p(p-1), but
    // still safer to use ints)
    for (int p=2; p<=(int)mN; p++)
    {
        for (int q=0; q<=(int)mN-p; q++)
        {
            ret += mAlpha[p][q] * p * (p-1) * pow(I1-3,(p-1)*(p-2)) * pow(I2-3,q);
        }
    }
    return ret;
}

double PolynomialMaterialLaw3d::Get_d2W_dI2(double I1, double I2)
{
    double ret = 0.0;

    // notes: use ints not unsigned as doing q-1
    // (except indexing from q=2 because multiplying by q(q-1), but
    // still safer to use ints)
    for (int p=0; p<=(int)mN; p++)
    {
        for (int q=2; q<=(int)mN-p; q++)
        {
            ret += mAlpha[p][q] * q * (q-1) * pow(I1-3,p) * pow(I2-3,(q-1)*(q-2));
        }
    }
    return ret;
}

double PolynomialMaterialLaw3d::Get_d2W_dI1I2(double I1, double I2)
{
    double ret = 0.0;

    // notes: use ints not unsigned as doing p-1
    // (except indexing from p=1,q=1 because multiplying by pq, but
    // still safer to use ints)
    for (int p=1; p<=(int)mN; p++)
    {
        for (int q=1; q<=(int)mN-p; q++)
        {
            ret += mAlpha[p][q] * p * q * pow(I1-3,p-1) * pow(I2-3,q-1);
        }
    }
    return ret;
}

double PolynomialMaterialLaw3d::GetAlpha(unsigned i, unsigned j)
{
    assert(i+j > 0);
    assert(i+j <= mN);

    return mAlpha[i][j];
}

PolynomialMaterialLaw3d::PolynomialMaterialLaw3d(unsigned n, std::vector<std::vector<double> > alpha)
{
    if (n==0)
    {
        EXCEPTION("n must be positive");
    }

    mN = n;

    // error checking: must have alpha[p][q]=0 if p+q>n
    for (unsigned p=0; p<=mN; p++)
    {
        if (alpha[p].size() < mN+1-p)
        {
            EXCEPTION("alpha not big enough");
        }

        for (unsigned q=0; q<alpha[p].size(); q++)
        {
            if ((p+q>mN) && (fabs(alpha[p][q]) > 1e-12))
            {
                EXCEPTION("alpha[" << p << "][" << q << "] should be zero, as p+q > " << n);
            }
        }
    }

    mAlpha = alpha;
}

std::vector<std::vector<double> > PolynomialMaterialLaw3d::GetZeroedAlpha(unsigned n)
{
    std::vector<std::vector<double> > alpha(n+1);

    for (unsigned i=0; i<n+1; i++)
    {
        alpha[i].resize(n+1);
        for (unsigned j=0; j<n+1; j++)
        {
            alpha[i][j] = 0.0;
        }
    }

    return alpha;
}
