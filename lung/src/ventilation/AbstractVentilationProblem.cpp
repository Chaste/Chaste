/*

Copyright (c) 2005-2015, University of Oxford.
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

#include "AbstractVentilationProblem.hpp"
#include "TrianglesMeshReader.hpp"
#include "MathsCustomFunctions.hpp"

AbstractVentilationProblem::AbstractVentilationProblem(const std::string& rMeshDirFilePath, unsigned rootIndex)
    : mOutletNodeIndex(rootIndex),
      mViscosity(1.92e-5),
      mDensity(1.15),
      mLengthScaling(1.0),
      mDynamicResistance(false),
      mRadiusOnEdge(false)
{
    TrianglesMeshReader<1,3> mesh_reader(rMeshDirFilePath);
    mMesh.ConstructFromMeshReader(mesh_reader);

    if (mMesh.GetNode(mOutletNodeIndex)->IsBoundaryNode() == false)
    {
        EXCEPTION("Outlet node is not a boundary node");
    }
}

void
AbstractVentilationProblem::SetMeshInMilliMetres()
{
    mLengthScaling = 1e-3;
}

void AbstractVentilationProblem::SetRadiusOnEdge(bool isOnEdges)
{
    mRadiusOnEdge = isOnEdges;
}

TetrahedralMesh<1, 3>&
AbstractVentilationProblem::rGetMesh()
{
    return mMesh;
}

double AbstractVentilationProblem::CalculateResistance(Element<1,3>& rElement, bool usePedley, double flux)
{
    //Resistance is based on radius, length and viscosity
    double radius = 0.0;
    if (mRadiusOnEdge)
    {
        radius = rElement.GetAttribute();
    }
    else
    {
        radius = ( rElement.GetNode(0)->rGetNodeAttributes()[0] + rElement.GetNode(1)->rGetNodeAttributes()[0]) / 2.0;
    }

    c_vector<double, 3> dummy;
    double length;
    mMesh.GetWeightedDirectionForElement(rElement.GetIndex(), dummy, length);

    radius *= mLengthScaling;
    length *= mLengthScaling;

    double resistance = 8.0*mViscosity*length/(M_PI*SmallPow(radius, 4));
    if ( usePedley )
    {
        /* Pedley et al. 1970
         * http://dx.doi.org/10.1016/0034-5687(70)90094-0
         * also Swan et al. 2012. 10.1016/j.jtbi.2012.01.042 (page 224)
         * Standard Poiseuille equation is similar to Pedley's modified Eq1. and matches Swan Eq5.
         * R_p = 128*mu*L/(pi*d^4) = 8*mu*L/(pi*r^4)
         *
         * Pedley Eq 2 and Swan Eq4:
         *  Z = C/(4*sqrt(2)) * sqrt(Re*d/l) = (C/4)*sqrt(Re*r/l)
         * Pedley suggests that C = 1.85
         * R_r = Z*R_p
         *
         * Reynold's number in a pipe is
         * Re = rho*v*d/mu (where d is a characteristic length scale - diameter of pipe)
         * since flux = v*area
         * Re = Q * d/(mu*area) = 2*rho*Q/(mu*pi*r) ... (see Swan p 224)
         *
         *
         * The upshot of this calculation is that the resistance is scaled with sqrt(Q)
         */
        double reynolds_number = fabs( 2.0 * mDensity * flux / (mViscosity * M_PI * radius) );
        double c = 1.85;
        double z = (c/4.0) * sqrt(reynolds_number * radius / length);
        // Pedley's method will only increase the resistance
        if (z > 1.0)
        {
            resistance *= z;
        }
    }
    return resistance;
}
