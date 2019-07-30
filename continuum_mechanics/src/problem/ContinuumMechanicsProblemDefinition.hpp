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

#ifndef CONTINUUMMECHANICSPROBLEMDEFINITION_HPP_
#define CONTINUUMMECHANICSPROBLEMDEFINITION_HPP_

#include "AbstractTetrahedralMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "AbstractIncompressibleMaterialLaw.hpp"
#include "AbstractCompressibleMaterialLaw.hpp"
#include "CompressibilityType.hpp"


/**
 *  Simple enumeration for denoting the type of body force, constant or defined via a function
 */
typedef enum BodyForceType_
{
    CONSTANT_BODY_FORCE,
    FUNCTIONAL_BODY_FORCE
} BodyForceType;

/**
 *  Simple enumeration for denoting the type of traction (Neumann) boundary condition:
 */
typedef enum TractionBoundaryConditionType_
{
    NO_TRACTIONS,
    ELEMENTWISE_TRACTION,
    FUNCTIONAL_TRACTION,
    PRESSURE_ON_DEFORMED,
    FUNCTIONAL_PRESSURE_ON_DEFORMED
} TractionBoundaryConditionType;


/**
 *  A class for specifying various parts of a continuum mechanics problem, Dirichlet node information (which
 *  nodes are in space in a solid problem, which nodes have fixed flow in a fluids problem), the body force (per unit mass)
 *  (usually acceleration due to gravity or zero), the traction boundary conditions, and the density.
 */
template<unsigned DIM>
class ContinuumMechanicsProblemDefinition
{
public:
    /** Special value for Dirichlet nodes, indicating that a Dirichlet boundary condition
     *  in a particular dimension is not specified */
    static const double FREE; // set to DBL_MAX

protected:
    /** The mesh being solved on */
    AbstractTetrahedralMesh<DIM,DIM>& mrMesh;

    /** Density of the body (constant throughout body) */
    double mDensity;

    //////////////////////////////
    // body force
    //////////////////////////////

    /** The body force type */
    BodyForceType mBodyForceType;

    /** The constant body force, only used if mBodyForceType is set appropriately */
    c_vector<double,DIM> mConstantBodyForce;

    /** The body force as a function of space and time, only used if mBodyForceType is set appropriately */
    c_vector<double,DIM> (*mpBodyForceFunction)(c_vector<double,DIM>& rX, double t);

    //////////////////////////////
    // tractions
    //////////////////////////////

    /** The traction (Neumann) boundary condition type */
    TractionBoundaryConditionType mTractionBoundaryConditionType;

    /** The surface elements on which tractions are applied */
    std::vector<BoundaryElement<DIM-1,DIM>*> mTractionBoundaryElements;

    /** The tractions on each surface element (only used if mTractionBoundaryConditionType is set appropriately) */
    std::vector<c_vector<double,DIM> > mElementwiseTractions;

    /** If the tractions are specified to correspond to a pressure acting on the surface: the pressure for the given
     *  boundary elements (only used if mTractionBoundaryConditionType is set appropriately) */
    double mNormalPressure;

    /** If the user asks this class to increment the pressure, the variable mNormalPressure will be altered
     *  depending on which increment it is. Here we store the original (full) pressure.
     */
    double mOriginalNormalPressure;

    /** The tractions as a function of space and time (only used if mTractionBoundaryConditionType is set appropriately) */
    c_vector<double,DIM> (*mpTractionBoundaryConditionFunction)(c_vector<double,DIM>& rX, double t);

    /** The normal pressure as a function if time (only used if mTractionBoundaryConditionType is set appropriately) */
    double (*mpNormalPressureFunction)(double t);


    ///////////////////////////////////////////
    // Dirichlet boundary conditions
    ///////////////////////////////////////////

    /**
     * All nodes (including non-vertices) which have a dirichlet boundary condition (ie position
     * prescribed in solid mechanics problems, flow prescribed in fluids problems). */
    std::vector<unsigned> mDirichletNodes;

    /** The values at the nodes with Dirichlet boundary conditions (displacement  */
    std::vector<c_vector<double,DIM> > mDirichletNodeValues;

    /** Whether the solver will be verbose or not. See dox for Set method below */
    bool mVerboseDuringSolve;

public:
    /**
     * Constructor initialises the body force to zero and density to 1.0
     * @param rMesh  is the mesh being solved on
     */
    ContinuumMechanicsProblemDefinition(AbstractTetrahedralMesh<DIM,DIM>& rMesh);

    /** Destructor */
    virtual ~ContinuumMechanicsProblemDefinition()
    {
    }

    /**
     *  Set the density
     *  @param density
     */
    void SetDensity(double density);

    /**
     * @return the density
     */
    double GetDensity();

    /**
     * Set a list of nodes (indices) to be given zero Dirichlet boundary condition
     * @param rZeroDirichletNodes the nodes at which the value (displacement/flow) is zero
     */
    void SetZeroDirichletNodes(std::vector<unsigned>& rZeroDirichletNodes);

    // No method here for setting non-zero Dirichlet BCs - these are in the subclasses..

    /**
     *  @return the Dirichlet nodes
     */
    std::vector<unsigned>& rGetDirichletNodes();

    /**
     * @return the Dirichlet node values.
     */
    std::vector<c_vector<double,DIM> >& rGetDirichletNodeValues();


    /**
     * Set the body force to be used - this version sets a constant body force
     * @param bodyForce the constant body force
     */
    void SetBodyForce(c_vector<double,DIM> bodyForce);

    /**
     * Set the body force to be used - this version sets a functional body force
     * @param pFunction a function of space and time returning a vector
     */
    void SetBodyForce(c_vector<double,DIM> (*pFunction)(c_vector<double,DIM>& rX, double t));

    /**
     * @return the body force at a particular point and time.
     * Note: The user can either call this, or check what type of body force has been set using
     * GetBodyForceType() and then call GetConstantBodyForce() or EvaluateBodyForceFunction(X,t).
     *
     * @param rX spatial location
     * @param t current time
     */
    c_vector<double,DIM> GetBodyForce(c_vector<double,DIM>& rX, double t = 0.0);

    /**
     * @return the body force type
     */
    BodyForceType GetBodyForceType();

    /**
     * @return the constant body force (error if GetBodyForceType()!=CONSTANT_BODY_FORCE)
     */
    c_vector<double,DIM> GetConstantBodyForce();

    /**
     * @return the body force function (error if GetBodyForceType()!=FUNCTIONAL_BODY_FORCE)
     *
     * @param rX spatial location
     * @param t current time
     */
    c_vector<double,DIM> EvaluateBodyForceFunction(c_vector<double,DIM>& rX, double t);

    /**
     * @return the traction (Neumann) boundary condition type
     */
    TractionBoundaryConditionType GetTractionBoundaryConditionType();

    /**
     * Set traction (Neumann) boundary conditions. This version takes in a vector of
     *  boundary elements, and corresponding tractions for each one.
     * @param rTractionBoundaryElements the boundary elements
     * @param rElementwiseTractions corresponding tractions
     */
    void SetTractionBoundaryConditions(std::vector<BoundaryElement<DIM-1,DIM>*>& rTractionBoundaryElements,
                                       std::vector<c_vector<double,DIM> >& rElementwiseTractions);

    /**
     * Set traction (Neumann) boundary conditions. This version takes in a vector of
     *  boundary elements, and a function to be evaluated at points in these boundary
     *  elements
     * @param rTractionBoundaryElements the boundary elements
     * @param pFunction the traction function (a function of space and time, returning a vector)
     */
    void SetTractionBoundaryConditions(std::vector<BoundaryElement<DIM-1,DIM>*>& rTractionBoundaryElements,
                                       c_vector<double,DIM> (*pFunction)(c_vector<double,DIM>& rX, double t));

    /**
     * Set traction (Neumann) boundary conditions. This version says that pressures should be applied
     * on surfaces in the DEFORMED body in the outward normal direction.
     *
     * @param rTractionBoundaryElements The boundary elements
     * @param normalPressure the corresponding pressure
     */
    void SetApplyNormalPressureOnDeformedSurface(std::vector<BoundaryElement<DIM-1,DIM>*>& rTractionBoundaryElements,
                                                 double normalPressure);

    /**
     * Set traction (Neumann) boundary conditions. This version says that pressures should be applied
     * on surfaces in the DEFORMED body in the outward normal direction, and here the pressure is specified
     * in FUNCTIONAL FORM
     *
     * @param rTractionBoundaryElements The boundary elements
     * @param pFunction the pressure function (a function of time, returning a double)
     */
    void SetApplyNormalPressureOnDeformedSurface(std::vector<BoundaryElement<DIM-1,DIM>*>& rTractionBoundaryElements,
                                                 double (*pFunction)(double t));

    /**
     * @return the vector of traction boundary elements
     */
    std::vector<BoundaryElement<DIM-1,DIM>*>& rGetTractionBoundaryElements();

    /**
     *  @return the element-wise tractions vector (corresponding to
     *  vector returned by rGetTractionBoundaryElements())
     *  (error if GetTractionBoundaryConditionType()!=ELEMENTWISE_TRACTION)
     */
    std::vector<c_vector<double,DIM> >& rGetElementwiseTractions();

    /**
     *  @return the pressure for the boundary elements (corresponding to
     *  vector returned by rGetTractionBoundaryElements())
     *  (error if GetTractionBoundaryConditionType()!=PRESSURE_ON_DEFORMED)
     */
    double GetNormalPressure();

    /**
     * Set the value that will be returned by GetNormalPressure() to be a fraction
     * of its full value.
     *
     * Note: you don't have to take into account previous calls when calling this.
     *  SetPressureScaling(0.1);
     *  SetPressureScaling(1);
     * will lead to the original pressure being used
     *
     * (Error if GetTractionBoundaryConditionType()!=PRESSURE_ON_DEFORMED)
     *
     * @param scaleFactor scale factor
     */
    void SetPressureScaling(double scaleFactor);

    /**
     * @return the traction boundary condition function (error if GetTractionBoundaryConditionType()!=FUNCTIONAL_TRACTION)
     *
     * @param rX spatial location
     * @param t current time
     */
    c_vector<double,DIM> EvaluateTractionFunction(c_vector<double,DIM>& rX, double t);

    /**
     * @return the pressure boundary condition function (error if GetTractionBoundaryConditionType()!=FUNCTIONAL_PRESSURE_ON_DEFORMED)
     *
     * @param t current time
     */
    double EvaluateNormalPressureFunction(double t);


    /**
     * Check all variables are set appropriately. Exceptions are thrown if any are not.
     * Derived classes can override but should call this version as well.
     */
    virtual void Validate();

    /**
     * Tell the solver to be verbose (print details on how the solve is progressing), or not.
     * @param verboseDuringSolve be verbose or not.
     */
    void SetVerboseDuringSolve(bool verboseDuringSolve = true)
    {
        mVerboseDuringSolve = verboseDuringSolve;
    }

    /**
     *  @return whether the solver should be verbose or not
     */
    bool GetVerboseDuringSolve()
    {
        return mVerboseDuringSolve;
    }
};

#endif // CONTINUUMMECHANICSPROBLEMDEFINITION_HPP_
