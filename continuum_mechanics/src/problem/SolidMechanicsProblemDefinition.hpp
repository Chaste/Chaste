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


#ifndef SOLIDMECHANICSPROBLEMDEFINITION_HPP_
#define SOLIDMECHANICSPROBLEMDEFINITION_HPP_

#include "ContinuumMechanicsProblemDefinition.hpp"
#include "QuadraticMesh.hpp"

/**
 *  A class for specifying various parts of a solid mechanics problem, in particular the material
 *  laws for the deforming body, and (inheriting functionality from a base class): fixed nodes information,
 *  the body force (per unit mass) (usually acceleration due to gravity or zero), the traction
 *  boundary conditions, and the density.
 */
template<unsigned DIM>
class SolidMechanicsProblemDefinition : public ContinuumMechanicsProblemDefinition<DIM>
{
private:
    /////////////////////////////
    // material law
    /////////////////////////////
    /**
     *  The material law, in the case of incompressible material laws. This vector is either of size 1, representing a
     *  homogeneous material, or of size num_elements, representing a heterogeneous material, with a material law per element.
     *  If he material is compressible, this vector will be of size zero.
     */
    std::vector<AbstractIncompressibleMaterialLaw<DIM>*> mIncompressibleMaterialLaws;
    /**
     *  The material law, in the case of compressible material laws. This vector is either of size 1, representing a
     *  homogeneous material, or of size num_elements, representing a heterogeneous material, with a material law per element.
     *  If the material is incompressible, this vector will be of size zero.
     */
    std::vector<AbstractCompressibleMaterialLaw<DIM>*> mCompressibleMaterialLaws;

    /** Whether the material is homogeneous (same material law everywhere) or heterogeneous */
    bool mIsHomogeneousMaterial;

    /** Whether the material is incompressible or compressible. (CompressibilityType is an enumeration).  */
    CompressibilityType mCompressibilityType;

    /** Whether the solver will use Petsc SNES or not. See dox for Set method below */
    bool mSolveUsingSnes;

    /**
     *  Helper function for checking whether a dynamic_cast succeeded or not, and throwing an exception
     *  if it failed.
     *
     *  @param compressibilityType compressibility type
     *  @param pMaterialLaw material law
     */
    void CheckCastSuccess(CompressibilityType compressibilityType,AbstractMaterialLaw<DIM>* pMaterialLaw);

public:
    /**
     * Constructor. Note body force initialised to zero and density to 1.0
     * @param rMesh Tesh being solved on
     */
    SolidMechanicsProblemDefinition(AbstractTetrahedralMesh<DIM,DIM>& rMesh);

    /** Destructor */
    virtual ~SolidMechanicsProblemDefinition()
    {
    }


    /**
     * Set a material law for the entire body (ie the homogeneous case). If compressibilityType==INCOMPRESSIBLE,
     * the material law pointer will be checked at run-time that it is of type `AbstractIncompressibleMaterialLaw`,
     * and similarly for the compressible case. Any previous material information will be deleted.
     *
     * @param compressibilityType either 'INCOMPRESSIBLE' or 'COMPRESSIBLE'
     * @param pMaterialLaw The material law for the entire body
     */
    void SetMaterialLaw(CompressibilityType compressibilityType, AbstractMaterialLaw<DIM>* pMaterialLaw);

    /**
     * Set a vector of material laws for the body, one for each element in the mesh (the heterogeneous case). If
     * compressibilityType==INCOMPRESSIBLE, the material law pointer will be checked at run-time that it is
     * of type `AbstractIncompressibleMaterialLaw`, and similarly for the compressible case.
     * Any previous material information will be deleted.
     *
     * @param compressibilityType either 'INCOMPRESSIBLE' or 'COMPRESSIBLE'
     * @param rMaterialLaws Vector of pointers to material laws
     */
    void SetMaterialLaw(CompressibilityType compressibilityType, std::vector<AbstractMaterialLaw<DIM>*>& rMaterialLaws);

    /**
     * @return whether the material is homogeneous or heterogeneous.
     * SetMaterialLaw() must be called before calling this.
     */
    bool IsHomogeneousMaterial();


    /**
     * @return whether the material is incompressible or compressible.
     * SetMaterialLaw() must be called before calling this. (Which can be checked
     * by calling Validate()).
     */
    CompressibilityType GetCompressibilityType();

    /**
     * @return the material law for a given element, when the body is incompressible. An assertion will
     * fail if GetCompressibilityType()!=INCOMPRESSIBLE. If the material is homogeneous, it doesn't matter what
     * the element index is.
     * @param elementIndex index of element
     */
    AbstractIncompressibleMaterialLaw<DIM>* GetIncompressibleMaterialLaw(unsigned elementIndex);

    /**
     * @return the material law for a given element, when the body is compressible. An assertion will
     * fail if GetCompressibilityType()!=COMPRESSIBLE. If the material is homogeneous, it doesn't matter what
     * the element index is.
     * @param elementIndex index of element
     */
    AbstractCompressibleMaterialLaw<DIM>* GetCompressibleMaterialLaw(unsigned elementIndex);

    /**
     * Set a list of nodes (indices) to be fixed in space with zero displacement
     * @param rFixedNodes the fixed nodes
     */
    void SetZeroDisplacementNodes(std::vector<unsigned>& rFixedNodes)
    {
        this->SetZeroDirichletNodes(rFixedNodes);
    }

    /**
     * Set a list of nodes to be fixed, with their corresponding new LOCATIONS (not displacements).
     * (This class will store as displacements though, and it is displacements that will be returned
     * by rGetDirichletNodeValues).
     * @param rFixedNodes the fixed node indices
     * @param rFixedNodeLocation corresponding locations
     */
    void SetFixedNodes(std::vector<unsigned>& rFixedNodes, std::vector<c_vector<double,DIM> >& rFixedNodeLocation);

    /**
     * Check all variables are set appropriately. Exceptions are thrown if any are not.
     * Derived classes can override but should call this version as well.
     */
    virtual void Validate();

    ///////////////////////////////////////////////////////////////////////////
    // The following methods set parameters used by the solver class
    // (AbstractNonlinearElasticitySolver). It is not ideal that they are
    // stored in this class - it would be nicer if this class was just about
    // the problem, not how it is solver - but in the absence of a globally
    // accessible config class, this is the best place.
    ///////////////////////////////////////////////////////////////////////////

    /**
     * Tell the solver class whether to use the PETSc SNES solver (the petsc nonlinear solver) or its
     * own nonlinear solve implementation.
     * @param solveUsingSnes solve using Snes or not
     */
    void SetSolveUsingSnes(bool solveUsingSnes = true)
    {
        mSolveUsingSnes = solveUsingSnes;
    }

    /**
     *  @return whether solver should use PETSc SNES nonlinear solver or not
     */
    bool GetSolveUsingSnes()
    {
        return mSolveUsingSnes;
    }
};

#endif /* SOLIDMECHANICSPROBLEMDEFINITION_HPP_ */
