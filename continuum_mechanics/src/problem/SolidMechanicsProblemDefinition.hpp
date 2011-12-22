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


#ifndef SOLIDMECHANICSPROBLEMDEFINITION_HPP_
#define SOLIDMECHANICSPROBLEMDEFINITION_HPP_

#include "ContinuumMechanicsProblemDefinition.hpp"

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
    SolidMechanicsProblemDefinition(QuadraticMesh<DIM>& rMesh);

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
     * Get whether the material is homogeneous or heterogeneous.
     * SetMaterialLaw() must be called before calling this.
     */
    bool IsHomogeneousMaterial();


    /**
     * Get whether the material is incompressible or compressible.
     * SetMaterialLaw() must be called before calling this. (Which can be checked
     * by calling Validate()).
     */
    CompressibilityType GetCompressibilityType();

    /**
     * Get the material law for a given element, when the body is incompressible. An assertion will
     * fail if GetCompressibilityType()!=INCOMPRESSIBLE. If the material is homogeneous, it doesn't matter what
     * the element index is.
     * @param elementIndex index of element
     */
    AbstractIncompressibleMaterialLaw<DIM>* GetIncompressibleMaterialLaw(unsigned elementIndex);

    /**
     * Get the material law for a given element, when the body is compressible. An assertion will
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

};




#endif /* SOLIDMECHANICSPROBLEMDEFINITION_HPP_ */
