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

#include "SimpleImpedanceProblem.hpp"
#include "TrianglesMeshReader.hpp"
#include "ReplicatableVector.hpp"

#include <cmath>

SimpleImpedanceProblem::SimpleImpedanceProblem(TetrahedralMesh<1,3>& rAirwaysMesh, unsigned rootIndex)
    : mrMesh(rAirwaysMesh),
      mOutletNodeIndex(rootIndex),
      mWalker(rAirwaysMesh, rootIndex),
      mRho(1.1500),                   //air density in Kg/m^3
      mMu(1.9e-5),                    //air viscosity in Pa s
      mH(5.8*98.0665*1e3),            //Tissue elastance in Pa/m^3 (5.8 cmH2O/L)
      mLengthScaling(1.0)
{
    mAcinarH = mH*(mrMesh.GetNumBoundaryNodes() - 1);

    mFrequencies.push_back(1.0);
    mFrequencies.push_back(2.0);
    mFrequencies.push_back(3.0);
    mFrequencies.push_back(5.0);
    mFrequencies.push_back(10.0);
    mFrequencies.push_back(20.0);
    mFrequencies.push_back(30.0);
}

SimpleImpedanceProblem::~SimpleImpedanceProblem()
{

}

void SimpleImpedanceProblem::SetFrequency(double frequency)
{
    mFrequencies.resize(1);
    mFrequencies[0] = frequency;
}

void SimpleImpedanceProblem::SetFrequencies(std::vector<double> frequencies)
{
    mFrequencies = frequencies;
}

std::vector<double>& SimpleImpedanceProblem::rGetFrequencies()
{
    return mFrequencies;
}

void SimpleImpedanceProblem::SetElastance(double elastance)
{
    mH = elastance;
    mAcinarH = mH*(mrMesh.GetNumBoundaryNodes() - 1);
}

void SimpleImpedanceProblem::Solve()
{
    Node<3>* p_node = mrMesh.GetNode(mOutletNodeIndex);
    Element<1,3>* p_element = mrMesh.GetElement(*(p_node->ContainingElementsBegin()));

    mImpedances.resize(mFrequencies.size());

    for (unsigned frequency_index = 0; frequency_index < mFrequencies.size(); ++frequency_index)
    {
        mImpedances[frequency_index] = CalculateElementImpedance(p_element, mFrequencies[frequency_index]);
    }
}

std::complex<double> SimpleImpedanceProblem::GetImpedance()
{
    return mImpedances[0];
}

std::vector<std::complex<double> >& SimpleImpedanceProblem::rGetImpedances()
{
    return mImpedances;
}


TetrahedralMesh<1, 3>& SimpleImpedanceProblem::rGetMesh()
{
    return mrMesh;
}

std::complex<double> SimpleImpedanceProblem::CalculateElementImpedance(Element<1,3>* pElement, double frequency)
{
    std::complex<double> Z(0, 0);

    //Get children
    if (mWalker.GetNumberOfChildElements(pElement) == 0u) //Branch is terminal, hence consider to be an acinus
    {
        Z = CalculateAcinusImpedance(mWalker.GetDistalNode(pElement), frequency);
    }
    else
    {
        std::complex<double> sum_one_over_Z_child(0,0); //Add up impedance of child elements

        std::vector<Element<1,3>* > child_eles = mWalker.GetChildElements(pElement);

        assert(child_eles.size() <= 2u);

        for (unsigned i = 0; i < child_eles.size(); ++i)
        {
            assert(child_eles[i] != pElement);

            std::complex<double> ele_impedance = CalculateElementImpedance(child_eles[i], frequency);

            if (real(ele_impedance) != 0.0 || imag(ele_impedance) != 0.0)
            {
                sum_one_over_Z_child += 1.0/ele_impedance;
            }
        }

        if (real(sum_one_over_Z_child) != 0.0 || imag(sum_one_over_Z_child) != 0.0)
        {
            Z = 1.0/sum_one_over_Z_child;
        }
    }


    double radius = (pElement->GetNode(0)->rGetNodeAttributes()[0] + pElement->GetNode(1)->rGetNodeAttributes()[0])/2.0; //Use average radius
    radius *= mLengthScaling;

    //For a 1D in 3D mesh, the element determinant == the element length
    c_matrix<double, 3, 1> jacobian; //not used
    double length;
    pElement->CalculateJacobian(jacobian, length);
    length *= mLengthScaling;

    double R = CalculateElementResistance(radius, length);
    double I = CalculateElementInertance(radius, length);

    double omega = 2*M_PI*frequency;
    std::complex<double> I_inertance(0, omega*I);

    return R + I_inertance + Z;
}

std::complex<double> SimpleImpedanceProblem::CalculateAcinusImpedance(Node<3>* pNode, double frequency)
{
    if (frequency == 0.0)
    {
        return 0.0;
    }

    std::complex<double> i(0, 1);
    double omega = 2*M_PI*frequency;

    return -i*mAcinarH/omega;
}


void SimpleImpedanceProblem::SetRho(double rho)
{
    mRho = rho;
}

void  SimpleImpedanceProblem::SetMu(double mu)
{
    mMu = mu;
}

double SimpleImpedanceProblem::CalculateElementResistance(double radius, double length)
{
    return 8*mMu*length/(M_PI*radius*radius*radius*radius);
}

double SimpleImpedanceProblem::CalculateElementInertance(double radius, double length)
{
    return mRho*length/(M_PI*radius*radius);
}

void SimpleImpedanceProblem::SetMeshInMilliMetres()
{
    mLengthScaling = 1e-3;
}
