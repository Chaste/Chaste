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
    StokesFlowProblemDefinition(QuadraticMesh<DIM>& rMesh)
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
     * Get the viscosity. Exception thrown if this hasn't been set yet.
     */
    double GetViscosity()
    {
        if(mMu < 0.0)
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
