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

#ifndef STOKESFLOWPROBLEMDEFINITION_HPP_
#define STOKESFLOWPROBLEMDEFINITION_HPP_

#include "ContinuumMechanicsProblemDefinition.hpp"

/**
 *  Class for defining everything needed for a solving the Stokes' Flow equations.
 *  This class mainly allows for the setting of the viscosity. It inherits from
 *  a class which allows for the setting of body force, density and various
 *  types of boundary condition.
 */
template<unsigned DIM>
class StokesFlowProblemDefinition : public ContinuumMechanicsProblemDefinition<DIM>
{
private:
    /** Dynamic viscosity */
    double mMu;

public:
    /**
     * Constructor (initialises viscosity to -1 so can check if it is unset
     * @param rMesh Quadratic mesh
     */
    StokesFlowProblemDefinition(AbstractTetrahedralMesh<DIM,DIM>& rMesh)
        : ContinuumMechanicsProblemDefinition<DIM>(rMesh),
          mMu(-1.0)
    {
    }

    /** Destructor */
    virtual ~StokesFlowProblemDefinition()
    {
    }

    /**
     * Set the viscosity
     * @param mu viscosity
     */
    void SetViscosity(double mu)
    {
        assert(mu > 0.0);
        mMu = mu;
    }

    /**
     * @return the viscosity. Exception thrown if this hasn't been set yet.
     */
    double GetViscosity()
    {
        if (mMu < 0.0)
        {
            EXCEPTION("Viscosity hasn't been set yet (for the Stokes' flow problem)");
        }
        return mMu;
    }

    /**
     * Set nodes on which to apply the boundary condition (u,v,w)=0, ie flow = 0 in all
     * components. Just calls SetZeroDirichletNodes() on the parent class.
     * @param rZeroFlowNodes Nodes on which to apply the boundary conditions
     */
    void SetZeroFlowNodes(std::vector<unsigned>& rZeroFlowNodes)
    {
        this->SetZeroDirichletNodes(rZeroFlowNodes);
    }

    /**
     * Set Dirichlet boundary conditions: provide a set of nodes and the prescribed flow
     * on these nodes.
     *
     * If you only want to set one component of the flow for the node rPrescribedFlowNodes[i],
     * set the other components of rPrescribedFlow[i] to be StokesFlowProblemDefinition::FREE.
     *
     * @param rPrescribedFlowNodes vector of node indices.
     * @param rPrescribedFlow vector of prescribed flow values for these nodes.
     */
    void SetPrescribedFlowNodes(std::vector<unsigned>& rPrescribedFlowNodes, std::vector<c_vector<double,DIM> >& rPrescribedFlow)
    {
        assert(rPrescribedFlowNodes.size()==rPrescribedFlow.size());

        this->mDirichletNodes = rPrescribedFlowNodes;
        this->mDirichletNodeValues = rPrescribedFlow;
    }

    /**
     * Check all variables are set appropriately. Exceptions are thrown if any are not.
     * Derived classes can override but should call this version as well.
     */
    virtual void Validate()
    {
        ContinuumMechanicsProblemDefinition<DIM>::Validate();
    }
};


#endif // STOKESFLOWPROBLEMDEFINITION_HPP_
