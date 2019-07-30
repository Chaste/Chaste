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

#ifndef MONODOMAINPURKINJECABLEASSEMBLER_HPP_
#define MONODOMAINPURKINJECABLEASSEMBLER_HPP_

#include "AbstractFeCableIntegralAssembler.hpp"
#include "HeartConfig.hpp"
#include "PdeSimulationTime.hpp"


/**
 * An assembler for the purkinje part of the left hand side matrix of the monodomain-purkinje linear system.
 * The full LHS matrix will look like (written in block form, although the code uses striping)
 * [ chi C M/dt + K   0                  ]
 * [      0           chi' C' M'/dt + K' ]
 * where the top left block is the standard LHS matrix in a monodomain problem, and the
 * bottom right block is the equivalent from integrals over cable elements.
 *
 * This class assembles the matrix
 * [  0          0           ]
 * [  0   chi' C' M'/dt + K' ]
 * The entries in this right-hand block are only non-zero on Purkinje nodes
 *
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MonodomainPurkinjeCableAssembler : public AbstractFeCableIntegralAssembler<ELEMENT_DIM,SPACE_DIM,2,false,true,NORMAL>
{
private:
    /**
     * Myocardium and purkinje potentials exist for each node, so PROBLEM_DIM=2.
     */
    static const unsigned PROBLEM_DIM=2;

    /**
     * Compute the cable element contribution to the matrix
     *
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases.
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i).
     * @param rX The point in space.
     * @param rU The unknown as a vector, u(i) = u_i.
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j).
     * @param pElement Pointer to the element.
     * @return the stencil matrix
     */
    c_matrix<double,PROBLEM_DIM*2,PROBLEM_DIM*2 /*2=number of bases per cable*/> ComputeCableMatrixTerm(
        c_vector<double, 2>& rPhi,
        c_matrix<double, ELEMENT_DIM, 2>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double,PROBLEM_DIM>& rU,
        c_matrix<double,PROBLEM_DIM, SPACE_DIM>& rGradU,
        Element<1,SPACE_DIM>* pElement)
    {
        c_matrix<double,PROBLEM_DIM*2, PROBLEM_DIM*2> ret = zero_matrix<double>(PROBLEM_DIM*2, PROBLEM_DIM*2);
        double capacitance = HeartConfig::Instance()->GetPurkinjeCapacitance();
        double chi = HeartConfig::Instance()->GetPurkinjeSurfaceAreaToVolumeRatio();
        double conductivity = HeartConfig::Instance()->GetPurkinjeConductivity();

        //We have to scale the assembled matrix by the cross sectional area of the Purkinje
        //fibre to ensure conservation of current at branch points. See #1899.
        const double fibre_cross_section_area = M_PI*pElement->GetAttribute()*pElement->GetAttribute();

        for(unsigned i=0; i<2; i++) // 2 = number of basis functions per cable element
        {
            for(unsigned j=0; j<2; j++)  // 2 = number of basis functions per cable element
            {
                ret(2*i,  2*j)   = 0;  // [V,V] block
                ret(2*i+1,2*j)   = 0;  // [Vpurkinje,V] block
                ret(2*i,  2*j+1) = 0;  // [V,Vpurkinje] block
                ret(2*i+1,2*j+1) = capacitance*chi*PdeSimulationTime::GetPdeTimeStepInverse()*rPhi(i)*rPhi(j);

                for (unsigned dim=0; dim<SPACE_DIM; dim++)
                {
                    ret(2*i+1,2*j+1) += conductivity*rGradPhi(dim,i)*rGradPhi(dim,j);
                }

                ret(2*i+1,2*j+1) *= fibre_cross_section_area;
            }
        }

        return ret;
    }
public:
    /**
     * Constructor
     * @param pMesh a pointer to a MixedDimensionMesh
     */
    MonodomainPurkinjeCableAssembler(MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
        : AbstractFeCableIntegralAssembler<ELEMENT_DIM,SPACE_DIM,2,false,true,NORMAL>(pMesh)
    {
        // Check radii have been set on the purkinje elements
        for (typename MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>::CableElementIterator iter = pMesh->GetCableElementIteratorBegin();
             iter != pMesh->GetCableElementIteratorEnd();
             ++iter)
        {
            double radius = (*iter)->GetAttribute();
            if (fabs(radius)<=DBL_EPSILON)
            {
                EXCEPTION("Radii not provided for all Purkinje elements - should be present in the mesh file or defined in test");
            }
        }
    }
};

#endif // MONODOMAINPURKINJECABLEASSEMBLER_HPP_
