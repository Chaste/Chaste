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

#ifndef ABSTRACTCONTINUUMMECHANICSSOLVER_HPP_
#define ABSTRACTCONTINUUMMECHANICSSOLVER_HPP_

#include "ContinuumMechanicsProblemDefinition.hpp"
#include "CompressibilityType.hpp"
#include "LinearSystem.hpp"
#include "OutputFileHandler.hpp"
#include "ReplicatableVector.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "QuadraticMesh.hpp"
#include "DistributedQuadraticMesh.hpp"
#include "Warnings.hpp"
#include "PetscException.hpp"
#include "GaussianQuadratureRule.hpp"
#include "PetscTools.hpp"
#include "MechanicsEventHandler.hpp"
#include "CommandLineArguments.hpp"
#include "VtkMeshWriter.hpp"


/**
 *  Simple enumeration for options that can be passed into
 *  AbstractContinuumMechanicsSolver::ApplyDirichletBoundaryConditions().
 *  See documentation for this method.
 */
typedef enum _ApplyDirichletBcsType
{
    LINEAR_PROBLEM,
    NONLINEAR_PROBLEM_APPLY_TO_RESIDUAL_ONLY,
    NONLINEAR_PROBLEM_APPLY_TO_EVERYTHING
} ApplyDirichletBcsType;

//forward declarations
template<unsigned DIM> class StressRecoveror;
template<unsigned DIM> class VtkNonlinearElasticitySolutionWriter;

/**
 *  General base class for continuum mechanics solvers. Deals with memory allocation,
 *  writing to file, applying boundary conditions, and adding an identity block (see below)
 *
 *  Note: the PROBLEM DIMENSION depends on DIM and whether a compressible or incompressible
 *  problem is being solved.
 *  (i) Compressible: mProblemDimension is set to DIM, and the ordering of the unknowns is, for
 *  2D say: [U1 V1 U2 V2 .. Un Vn]
 *  (ii) Incompressible: Here we have spatial unknowns (at all nodes of the QuadraticMesh) and
 *  pressure unknowns (only at the vertices, as linear bases are). We choose mProblemDim = DIM+1,
 *  use (for parallelisation reasons) the ordering [U1 V1 P1 , .. , Un Vn, Pn], introducing dummy
 *  variables for pressure unknowns at internal hences, with the equation Pi=0 at these nodes,
 *  hence the need for an identity block to be added the corresponding parts of the matrix.
 */
template<unsigned DIM>
class AbstractContinuumMechanicsSolver
{
    friend class StressRecoveror<DIM>;
    friend class VtkNonlinearElasticitySolutionWriter<DIM>;

protected:
    /**
     * The mesh to be solved on. Requires 6 nodes per triangle (or 10 per tetrahedron)
     * as quadratic bases are used.
     */
    AbstractTetrahedralMesh<DIM, DIM>& mrQuadMesh;

    /** Problem definition class - contains info on boundary conditions, etc */
    ContinuumMechanicsProblemDefinition<DIM>& mrProblemDefinition;

    /** Whether to write any output. */
    bool mWriteOutput;

    /** Where to write output, relative to CHASTE_TEST_OUTPUT. */
    std::string mOutputDirectory;

    /** Output file handler. */
    OutputFileHandler* mpOutputFileHandler;

    /** Spatial solution:
     *  For solids problems, mSpatialSolution[i](j) = x_j (new position) for node i.
     *  For fluids problems, mSpatialSolution[i](j) = u_j (flow) for node i. */
    std::vector<c_vector<double,DIM> > mSpatialSolution;

    /** Pressures solution at each vertex of the mesh. Only valid if mCompressibilityType==INCOMPRESSIBLE. */
    std::vector<double> mPressureSolution;


    /**
     * The current solution, in the form (assuming 2d):
     *   Incompressible problem: [u1 v1 u2 v2 ... uN vN p1 p2 .. pM]
     *   Compressible problem:   [u1 v1 u2 v2 ... uN vN]
     * where there are N total nodes and M vertices.
     */
    std::vector<double> mCurrentSolution;

    /** Gaussian quadrature rule. */
    GaussianQuadratureRule<DIM>* mpQuadratureRule;

    /** Boundary Gaussian quadrature rule. */
    GaussianQuadratureRule<DIM-1>* mpBoundaryQuadratureRule;

    /**
     * This is equal to either COMPRESSIBLE or INCOMPRESSIBLE (see enumeration class)
     * and is only used in computing mNumDofs and allocating matrix memory.
     */
    CompressibilityType mCompressibilityType;

    /**
     *  The problem dimension - the number of unknowns at each node.
     *  For compressible problems, where only the deformation/flow is computed at each node, this is
     *  equal to DIM
     *  For incompressible problems, where pressure is computed as well, this is equal to DIM+1
     *  (there is a pressure variable defined even for internal nodes at which pressure is not computed,
     *  this is a dummy pressure variable -- done like this for parallelisation reasons. See above).
     */
    unsigned mProblemDimension;

    /**
     * Number of degrees of freedom (equal to mProblemDim*num_nodes)
     */
    unsigned mNumDofs;

    /**
     * If SetVerboseDuringSolve() is called on the problem definition class, or the command line argument
     * "-mech_verbose" or "-mech_very_verbose" is given, than this bool will be  set to true and lots
     * of details about each nonlinear solve (including timing breakdowns) will be printed out
     */
    bool mVerbose;

    /**
     * Residual vector nonlinear problems.
     *
     * Since the residual in nonlinear problems is usually also the RHS vector in the linear
     * system, it may seem unncessary to also have the member variable mLinearSystemRhsVector.
     *
     * However: Newton's method is Ju = f, where J is the Jacobian, u the (negative of the) update
     * and f the residual, but when applying Dirichlet boundary conditions in
     * the compressible case, we alter the rows of the matrix and also alter the columns in order to
     * maintain symmetry. This requires making further changes to the right-hand vector, meaning that
     * it no longer properly represents the residual. Hence, we have to use two vectors.
     *
     * Overall, this can be represents as
     *  - compute residual f
     *  - compute Jacobian J
     *  - apply BCs to f.
     *  - alter the linear system from Ju=f to (J*)u=f* which enforces the dirichlet boundary conditions but enforces them symmetrically.
     *
     * mLinearSystemRhsVector below represents f*.
     */
    Vec mResidualVector;

    /**
     * The RHS side in the linear system that is solved each Newton iteration.
     */
    Vec mLinearSystemRhsVector;

    /**
     * Jacobian matrix of the nonlinear system, LHS matrix for the linear system.
     */
    Mat mSystemLhsMatrix;

    /**
     * Helper vector (see ApplyDirichletBoundaryConditions code).
     */
    Vec mDirichletBoundaryConditionsVector;

    /**
     * Precondition matrix for the linear system.
     */
    Mat mPreconditionMatrix;

    /**
     * Allocates memory for the matrices and vectors
     */
    void AllocateMatrixMemory();


    /**
     * Apply the Dirichlet boundary conditions to the linear system.
     *
     * The first input parameter should be one of the following
     *   LINEAR_PROBLEM -- indicating the overall problem is linear, and in which
     *     case the BCs will be applied to both the matrix and vector of the
     *     linear system
     *
     *   NONLINEAR_PROBLEM_APPLY_TO_RESIDUAL_ONLY -- indicating the overall problem is nonlinear
     *     and here only the residual vector will be altered (apply the Dirichlet boundary
     *     conditions to the residual vector involves setting appropriate components to the
     *     difference between the current value and the correct value).
     *
     *   NONLINEAR_PROBLEM_APPLY_TO_EVERYTHING -- indicating the overall problem is nonlinear,
     *     and here, the residual vector will be altered, as will the matrix and RHS vector
     *     (see documentation on mResidualVector for why there is a separate residual vector
     *     and RHS vector).
     *
     *  The second parameter should be true if the overall problem is symmetric, in which case
     *  boundary conditions will be applied keeping the matrix symmetric (if the matrix is being
     *  altered). See in-code comments for how this is done.
     *
     *  @param type see above
     *  @param symmetricProblem see above
     */
    void ApplyDirichletBoundaryConditions(ApplyDirichletBcsType type, bool symmetricProblem);

    /**
     * For incompressible problems, we use the following ordering:
     * [U0 V0 W0 P0 U1 V1 W1 P1 .. Un Vn Wn Pn]
     * where (U V W) is the displacement, and P is the pressure.
     * Therefore P is defined at every node; however for the solve, linear basis functions are used for
     * the pressure, so pressure is solved for at each vertex, not at internal nodes (which are for
     * quadratic basis functions).
     *
     * The pressure variables for non-vertex nodes are therefore dummy variables. This method enforces the
     * condition P_i=0, where i corresponds to a non-vertex node.
     *
     * The first input parameter should be one of the following
     *   LINEAR_PROBLEM -- indicating the overall problem is linear, and in which
     *     case the matrix will be altered (1 on diagonal, zeros assumed on rest of row)
     *     and the rhs vector will be set to 0.
     *
     *   NONLINEAR_PROBLEM_APPLY_TO_RESIDUAL_ONLY -- indicating the overall problem is nonlinear
     *     and here only the residual vector will be altered (sets the residual value of the appropriate
     *     rows to P_i-0.0.
     *
     *   NONLINEAR_PROBLEM_APPLY_TO_EVERYTHING -- indicating the overall problem is nonlinear,
     *     and here, the residual vector will be altered, as will the matrix and RHS vector
     *     (see documentation on mResidualVector for why there is a separate residual vector
     *     and RHS vector).
     *
     *  Note the row (and columns) for the dummy variables are not explicitly zeroed, as they
     *  are assumed to be zero already - in a finite element assembly nothing will have been
     *  written to the rows or columns corresponding to dummy variables. Hence this method
     *  always maintains symmetry.
     *
     *  @param type see above
     */
    void AddIdentityBlockForDummyPressureVariables(ApplyDirichletBcsType type);


    /**
     * For incompressible problems, we use the following ordering:
     * [U0 V0 W0 P0 U1 V1 W1 P1 .. Un Vn Wn Pn]
     * where (U V W) is the displacement, and P is the pressure.
     * Therefore P is defined at every node; however for the solve, linear basis functions are used for
     * the pressure, so pressure is solved for at each vertex, not at internal nodes (which are for
     * quadratic basis functions).
     *
     * The above method, AddIdentityBlockForDummyPressureVariables(), enforces the
     * condition P_i=0, where i corresponds to a non-vertex node.
     *
     * AFTER the solve, this method can be used to remove the dummy values, but looping over all
     * edges, and linearly interpolating the pressure at the two vertices onto the internal node.
     *
     * This method assumes each internal node is midway between the two vertices.
     */
    void RemovePressureDummyValuesThroughLinearInterpolation();

public:
    /**
     *  Constructor
     *  @param rQuadMesh the mesh
     *  @param rProblemDefinition problem definition object
     *  @param outputDirectory output directory name
     *  @param compressibilityType 'INCOMPRESSIBLE' or 'COMPRESSIBLE'
     */
    AbstractContinuumMechanicsSolver(AbstractTetrahedralMesh<DIM, DIM>& rQuadMesh,
                                     ContinuumMechanicsProblemDefinition<DIM>& rProblemDefinition,
                                     std::string outputDirectory,
                                     CompressibilityType compressibilityType);

    /** Destructor */
    virtual ~AbstractContinuumMechanicsSolver();

    /**
     * Write the spatial solution (deformed position if solids, flow if fluids) at the nodes
     *
     * @param fileName (stem)
     * @param fileExtension to append at end.
     * @param counterToAppend append a counter to the file name (defaults to nothing appended).
     *
     * For example:
     * WriteCurrentSpatialSolution("solution","nodes") --> file called "solution.nodes"
     * WriteCurrentSpatialSolution("solution","nodes",3) --> file called "solution_3.nodes"
     */
    void WriteCurrentSpatialSolution(std::string fileName, std::string fileExtension, int counterToAppend=-1);

    /**
     * Write the pressure solution. Only valid if mCompressibilityType==INCOMPRESSIBLE.
     * Writes the pressure for ALL nodes on the mesh, including internal nodes (as these are not assumed to
     * be have indices greater than vertex nodes). As linear basis functions are used for pressure, the
     * pressure solution is only computed at the vertices, and pressure dummy variables are used at internal
     * nodes, hence for each internal node, 0 will be written to file.
     *
     * @param counterToAppend append a counter to the file name
     *
     * WriteCurrentPressureSolution() --> file called "pressure.txt"
     * WriteCurrentPressureSolution(3) --> file called "pressure_3.txt"
     */
    void WriteCurrentPressureSolution( int counterToAppend=-1);

    /**
     * Set whether to write any output.
     *
     * @param writeOutput (defaults to true)
     */
    void SetWriteOutput(bool writeOutput=true);

    /**
     * Convert the output to vtk format (placed in a folder called vtk in the output directory).
     * @param spatialSolutionName is used to identify the spatial solution as a velocity, displacement...
     *
     * ** TO BE DEPRECATED - see #2321 **
     *
     */
    void CreateVtkOutput(std::string spatialSolutionName="Spatial solution");

    /**
     * @return the current solution vector (advanced use only - for getting the deformed position use
     * rGetDeformedPosition()).
     */
    std::vector<double>& rGetCurrentSolution()
    {
        return mCurrentSolution;
    }

    /**
     * @return the spatial solution. For solids problems this will be the deformed position,
     * for fluids problems this will be the flow.
     */
    virtual std::vector<c_vector<double,DIM> >& rGetSpatialSolution()=0;

    /**
     *  @return the pressure, for each NODE in the mesh. If the node is an internal node of the quadratic mesh
     *  the pressure is not computed at this node (as linear basis functions are used for the pressure, so
     *  pressures unknowns are only present at vertices), so a dummy pressure value of 0 is returned in this
     *  vector.
     *
     *  Only valid if mCompressibilityType==INCOMPRESSIBLE
     */
    std::vector<double>& rGetPressures();
};

template<unsigned DIM>
AbstractContinuumMechanicsSolver<DIM>::AbstractContinuumMechanicsSolver(AbstractTetrahedralMesh<DIM, DIM>& rQuadMesh,
                                                                        ContinuumMechanicsProblemDefinition<DIM>& rProblemDefinition,
                                                                        std::string outputDirectory,
                                                                        CompressibilityType compressibilityType)
    : mrQuadMesh(rQuadMesh),
      mrProblemDefinition(rProblemDefinition),
      mOutputDirectory(outputDirectory),
      mpOutputFileHandler(nullptr),
      mpQuadratureRule(nullptr),
      mpBoundaryQuadratureRule(nullptr),
      mCompressibilityType(compressibilityType),
      mResidualVector(nullptr),
      mSystemLhsMatrix(nullptr),
      mPreconditionMatrix(nullptr)
{
    assert(DIM==2 || DIM==3);

    //Check that the mesh is Quadratic
    QuadraticMesh<DIM>* p_quad_mesh = dynamic_cast<QuadraticMesh<DIM>* >(&rQuadMesh);
    DistributedQuadraticMesh<DIM>* p_distributed_quad_mesh = dynamic_cast<DistributedQuadraticMesh<DIM>* >(&rQuadMesh);

    if ((p_quad_mesh == nullptr) && (p_distributed_quad_mesh == nullptr))
    {
        EXCEPTION("Continuum mechanics solvers require a quadratic mesh");
    }


    mVerbose = (mrProblemDefinition.GetVerboseDuringSolve() ||
                CommandLineArguments::Instance()->OptionExists("-mech_verbose") ||
                CommandLineArguments::Instance()->OptionExists("-mech_very_verbose") );

    mWriteOutput = (mOutputDirectory != "");
    if (mWriteOutput)
    {
        mpOutputFileHandler = new OutputFileHandler(mOutputDirectory);
    }

    // See dox for mProblemDimension
    mProblemDimension = mCompressibilityType==COMPRESSIBLE ? DIM : DIM+1;
    mNumDofs = mProblemDimension*mrQuadMesh.GetNumNodes();

    AllocateMatrixMemory();

    // In general the Jacobian for a mechanics problem is non-polynomial.
    // We therefore use the highest order integration rule available.
    mpQuadratureRule = new GaussianQuadratureRule<DIM>(3);
    // The boundary forcing terms (or tractions) are also non-polynomial in general.
    // Again, we use the highest order integration rule available.
    mpBoundaryQuadratureRule = new GaussianQuadratureRule<DIM-1>(3);

    mCurrentSolution.resize(mNumDofs, 0.0);
}

template<unsigned DIM>
AbstractContinuumMechanicsSolver<DIM>::~AbstractContinuumMechanicsSolver()
{
    if (mpOutputFileHandler)
    {
        delete mpOutputFileHandler;
    }

    if (mpQuadratureRule)
    {
        delete mpQuadratureRule;
        delete mpBoundaryQuadratureRule;
    }

    if (mResidualVector)
    {
        PetscTools::Destroy(mResidualVector);
        PetscTools::Destroy(mLinearSystemRhsVector);
        PetscTools::Destroy(mSystemLhsMatrix);
        PetscTools::Destroy(mPreconditionMatrix);
    }

    if (mDirichletBoundaryConditionsVector)
    {
        PetscTools::Destroy(mDirichletBoundaryConditionsVector);
    }
}

template<unsigned DIM>
void AbstractContinuumMechanicsSolver<DIM>::WriteCurrentSpatialSolution(std::string fileName,
                                                                        std::string fileExtension,
                                                                        int counterToAppend)
{
    // Only write output if the flag mWriteOutput has been set
    if (!mWriteOutput)
    {
        return;
    }

    if (PetscTools::AmMaster())
    {
        std::stringstream file_name;

        file_name << fileName;
        if (counterToAppend >= 0)
        {
            file_name << "_" << counterToAppend;
        }
        file_name << "." << fileExtension;

        out_stream p_file = mpOutputFileHandler->OpenOutputFile(file_name.str());

        std::vector<c_vector<double,DIM> >& r_spatial_solution = rGetSpatialSolution();
        for (unsigned i=0; i<r_spatial_solution.size(); i++)
        {
    //        for (unsigned j=0; j<DIM; j++)
    //        {
    //            *p_file << mrQuadMesh.GetNode(i)->rGetLocation()[j] << " ";
    //        }
            for (unsigned j=0; j<DIM; j++)
            {
                *p_file << r_spatial_solution[i](j) << " ";
            }
            *p_file << "\n";
        }
        p_file->close();
    }
    PetscTools::Barrier("WriteSpatial");
}

template<unsigned DIM>
void AbstractContinuumMechanicsSolver<DIM>::WriteCurrentPressureSolution(int counterToAppend)
{
    // Only write output if the flag mWriteOutput has been set
    if (mWriteOutput)
    {
        // Only the master writes
        if (PetscTools::AmMaster() && mWriteOutput)
        {
            std::stringstream file_name;

            file_name << "pressure";
            if (counterToAppend >= 0)
            {
                file_name << "_" << counterToAppend;
            }
            file_name << ".txt";

            out_stream p_file = mpOutputFileHandler->OpenOutputFile(file_name.str());

            std::vector<double> &r_pressure = rGetPressures();
            for (unsigned i = 0; i < r_pressure.size(); i++)
            {
                for (unsigned j = 0; j < DIM; j++)
                {
                    *p_file << mrQuadMesh.GetNode(i)->rGetLocation()[j] << " ";
                }

                *p_file << r_pressure[i] << "\n";
            }
            p_file->close();
        }
        PetscTools::Barrier("WritePressure");
    }
}

template<unsigned DIM>
void AbstractContinuumMechanicsSolver<DIM>::SetWriteOutput(bool writeOutput)
{
    if (writeOutput && (mOutputDirectory==""))
    {
        EXCEPTION("Can't write output if no output directory was given in constructor");
    }
    mWriteOutput = writeOutput;
}

// ** TO BE DEPRECATED - see #2321 **
template<unsigned DIM>
void AbstractContinuumMechanicsSolver<DIM>::CreateVtkOutput(std::string spatialSolutionName)
{
    if (this->mOutputDirectory=="")
    {
        EXCEPTION("No output directory was given so no output was written, cannot convert to VTK format");
    }
#ifdef CHASTE_VTK
    VtkMeshWriter<DIM, DIM> mesh_writer(this->mOutputDirectory + "/vtk", "solution", true);

    mesh_writer.AddPointData(spatialSolutionName, this->rGetSpatialSolution());

    if (mCompressibilityType==INCOMPRESSIBLE)
    {
        mesh_writer.AddPointData("Pressure", rGetPressures());
    }

    //Output the element attribute as cell data.
    std::vector<double> element_attribute;
    for (typename QuadraticMesh<DIM>::ElementIterator iter = this->mrQuadMesh.GetElementIteratorBegin();
        iter != this->mrQuadMesh.GetElementIteratorEnd();
        ++iter)
    {
        element_attribute.push_back(iter->GetAttribute());
    }
    mesh_writer.AddCellData("Attribute", element_attribute);

    mesh_writer.WriteFilesUsingMesh(this->mrQuadMesh);
#endif
}

template<unsigned DIM>
std::vector<double>& AbstractContinuumMechanicsSolver<DIM>::rGetPressures()
{
    assert(mProblemDimension==DIM+1);

    mPressureSolution.clear();
    mPressureSolution.resize(mrQuadMesh.GetNumNodes());

    for (unsigned i=0; i<mrQuadMesh.GetNumNodes(); i++)
    {
        mPressureSolution[i] = mCurrentSolution[mProblemDimension*i + DIM];
    }
    return mPressureSolution;
}

template<unsigned DIM>
void AbstractContinuumMechanicsSolver<DIM>::RemovePressureDummyValuesThroughLinearInterpolation()
{
    assert(mProblemDimension==DIM+1);

    // For quadratic triangles, node 3 is between nodes 1 and 2, node 4 is between 0 and 2, etc
    unsigned internal_nodes_2d[3] = {3,4,5};
    unsigned neighbouring_vertices_2d[3][2] = { {1,2}, {2,0}, {0,1} };

    // ordering for quadratic tetrahedra
    unsigned internal_nodes_3d[6] = {4,5,6,7,8,9};
    unsigned neighbouring_vertices_3d[6][2] = { {0,1}, {1,2}, {0,2}, {0,3}, {1,3}, {2,3} };

    unsigned num_internal_nodes_per_element = DIM==2 ? 3 : 6;

    // loop over elements, then loop over edges.
    for (typename AbstractTetrahedralMesh<DIM,DIM>::ElementIterator iter = mrQuadMesh.GetElementIteratorBegin();
         iter != mrQuadMesh.GetElementIteratorEnd();
         ++iter)
    {
        for (unsigned i=0; i<num_internal_nodes_per_element; i++)
        {
            unsigned global_index;
            double left_val;
            double right_val;

            if (DIM == 2)
            {
                global_index = iter->GetNodeGlobalIndex( internal_nodes_2d[i] );
                unsigned vertex_0_global_index =iter->GetNodeGlobalIndex( neighbouring_vertices_2d[i][0] );
                unsigned vertex_1_global_index =iter->GetNodeGlobalIndex( neighbouring_vertices_2d[i][1] );
                left_val = mCurrentSolution[mProblemDimension*vertex_0_global_index + DIM];
                right_val = mCurrentSolution[mProblemDimension*vertex_1_global_index + DIM];
            }
            else
            {
                global_index = iter->GetNodeGlobalIndex( internal_nodes_3d[i] );
                unsigned vertex_0_global_index =iter->GetNodeGlobalIndex( neighbouring_vertices_3d[i][0] );
                unsigned vertex_1_global_index =iter->GetNodeGlobalIndex( neighbouring_vertices_3d[i][1] );
                left_val = mCurrentSolution[mProblemDimension*vertex_0_global_index + DIM];
                right_val = mCurrentSolution[mProblemDimension*vertex_1_global_index + DIM];
            }

            // this line assumes the internal node is midway between the two vertices
            mCurrentSolution[mProblemDimension*global_index + DIM] =  0.5 * (left_val + right_val);
        }
    }
}

/*
 * This method applies the appropriate BCs to the residual and/or linear system
 *
 * For the latter, and second input parameter==true, the BCs are imposed in such a way as
 * to ensure that a symmetric linear system remains symmetric. For each row with boundary condition
 * applied, both the row and column are zero'd and the RHS vector modified to take into account the
 * zero'd column.
 *
 * Suppose we have a matrix
 * [a b c] [x] = [ b1 ]
 * [d e f] [y]   [ b2 ]
 * [g h i] [z]   [ b3 ]
 * and we want to apply the boundary condition x=v without losing symmetry if the matrix is
 * symmetric. We apply the boundary condition
 * [1 0 0] [x] = [ v  ]
 * [d e f] [y]   [ b2 ]
 * [g h i] [z]   [ b3 ]
 * and then zero the column as well, adding a term to the RHS to take account for the
 * zero-matrix components
 * [1 0 0] [x] = [ v  ] - v[ 0 ]
 * [0 e f] [y]   [ b2 ]    [ d ]
 * [0 h i] [z]   [ b3 ]    [ g ]
 * Note the last term is the first column of the matrix, with one component zeroed, and
 * multiplied by the boundary condition value. This last term is then stored in
 * mDirichletBoundaryConditionsVector, and in general is equal to:
 * SUM_{d=1..D} v_d a'_d
 * where v_d is the boundary value of boundary condition d (d an index into the matrix),
 * and a'_d is the dth-column of the matrix but with the d-th component zeroed, and where
 * there are D boundary conditions
 */
template<unsigned DIM>
void AbstractContinuumMechanicsSolver<DIM>::ApplyDirichletBoundaryConditions(ApplyDirichletBcsType type, bool symmetricProblem)
{
    std::vector<unsigned> rows;
    std::vector<double> values;

    // Whether to apply symmetrically, ie alter columns as well as rows (see comment above)
    bool applySymmetrically = (type!=NONLINEAR_PROBLEM_APPLY_TO_RESIDUAL_ONLY) && symmetricProblem;

    if (applySymmetrically)
    {
        if (mDirichletBoundaryConditionsVector == nullptr)
        {
            VecDuplicate(mResidualVector, &mDirichletBoundaryConditionsVector);
        }

        PetscVecTools::Zero(mDirichletBoundaryConditionsVector);
        PetscMatTools::Finalise(mSystemLhsMatrix);
    }

    ///////////////////////////////////////
    // collect the entries to be altered
    ///////////////////////////////////////

    for (unsigned i=0; i<mrProblemDefinition.rGetDirichletNodes().size(); i++)
    {
        unsigned node_index = mrProblemDefinition.rGetDirichletNodes()[i];

        for (unsigned j=0; j<DIM; j++)
        {
            double dirichlet_val = mrProblemDefinition.rGetDirichletNodeValues()[i](j);

            if (dirichlet_val != ContinuumMechanicsProblemDefinition<DIM>::FREE)
            {
                double val;
                unsigned dof_index = mProblemDimension*node_index+j;

                if (type == LINEAR_PROBLEM)
                {
                    val = dirichlet_val;
                }
                else
                {
                    // The boundary conditions on the NONLINEAR SYSTEM are x=boundary_values
                    // on the boundary nodes. However:
                    // The boundary conditions on the LINEAR SYSTEM at each Newton step (Ju=f,
                    // where J is the Jacobian, u the negative update vector and f is the residual) is:
                    // u=current_soln-boundary_values on the boundary nodes
                    val = mCurrentSolution[dof_index] - dirichlet_val;
                }
                rows.push_back(dof_index);
                values.push_back(val);
            }
        }
    }

    ///////////////////////////////////////
    // do the alterations
    ///////////////////////////////////////

    if (applySymmetrically)
    {
        // Modify the matrix columns
        for (unsigned i=0; i<rows.size(); i++)
        {
            unsigned col = rows[i];
            double minus_value = -values[i];

            // Get a vector which will store the column of the matrix (column d, where d is
            // the index of the row (and column) to be altered for the boundary condition.
            // Since the matrix is symmetric when get row number "col" and treat it as a column.
            // PETSc uses compressed row format and therefore getting rows is far more efficient
            // than getting columns.
            Vec matrix_col = PetscMatTools::GetMatrixRowDistributed(mSystemLhsMatrix,col);

            // Zero the correct entry of the column
            PetscVecTools::SetElement(matrix_col, col, 0.0);

            // Set up the RHS Dirichlet boundary conditions vector
            // E.g. for a boundary node at the zeroth node (x_0 = value), this is equal to
            //   -value*[0 a_21 a_31 .. a_N1]
            // and will be added to the RHS.
            PetscVecTools::AddScaledVector(mDirichletBoundaryConditionsVector, matrix_col, minus_value);
            PetscTools::Destroy(matrix_col);
        }
    }

    if (type!=NONLINEAR_PROBLEM_APPLY_TO_RESIDUAL_ONLY) // i.e. doing a whole linear system
    {
        // Now zero the appropriate rows and columns of the matrix. If the matrix is symmetric we apply the
        // boundary conditions in a way the symmetry isn't lost (rows and columns). If not only the row is
        // zeroed.
        if (applySymmetrically)
        {
            PetscMatTools::ZeroRowsAndColumnsWithValueOnDiagonal(mSystemLhsMatrix, rows, 1.0);
            PetscMatTools::ZeroRowsAndColumnsWithValueOnDiagonal(mPreconditionMatrix, rows, 1.0);

            // Apply the RHS boundary conditions modification if required.
            PetscVecTools::AddScaledVector(mLinearSystemRhsVector, mDirichletBoundaryConditionsVector, 1.0);
        }
        else
        {
            PetscMatTools::ZeroRowsWithValueOnDiagonal(mSystemLhsMatrix, rows, 1.0);
            PetscMatTools::ZeroRowsWithValueOnDiagonal(mPreconditionMatrix, rows, 1.0);
        }
    }

    if (type!=LINEAR_PROBLEM)
    {
        for (unsigned i=0; i<rows.size(); i++)
        {
            PetscVecTools::SetElement(mResidualVector, rows[i], values[i]);
        }
    }

    if (type!=NONLINEAR_PROBLEM_APPLY_TO_RESIDUAL_ONLY)
    {
        for (unsigned i=0; i<rows.size(); i++)
        {
            PetscVecTools::SetElement(mLinearSystemRhsVector, rows[i], values[i]);
        }
    }
}

template<unsigned DIM>
void AbstractContinuumMechanicsSolver<DIM>::AddIdentityBlockForDummyPressureVariables(ApplyDirichletBcsType type)
{
    assert(mCompressibilityType==INCOMPRESSIBLE);

    int lo, hi;
    VecGetOwnershipRange(mResidualVector, &lo, &hi);

    for (unsigned i=0; i<mrQuadMesh.GetNumNodes(); i++)
    {
        if (mrQuadMesh.GetNode(i)->IsInternal())
        {
            unsigned row = (DIM+1)*i + DIM; // DIM+1 is the problem dimension
            if (lo <= (int)row && (int)row < hi)
            {
                if (type!=LINEAR_PROBLEM)
                {
                    PetscVecTools::SetElement(mResidualVector, row, mCurrentSolution[row]-0.0);
                }
                if (type!=NONLINEAR_PROBLEM_APPLY_TO_RESIDUAL_ONLY) // ie doing a whole linear system
                {
                    double rhs_vector_val = type==LINEAR_PROBLEM ? 0.0 : mCurrentSolution[row]-0.0;
                    PetscVecTools::SetElement(mLinearSystemRhsVector, row, rhs_vector_val);
                    // This assumes the row is already zero, which is should be..
                    PetscMatTools::SetElement(mSystemLhsMatrix, row, row, 1.0);
                    PetscMatTools::SetElement(mPreconditionMatrix, row, row, 1.0);
                }
            }
        }
    }
}

template<unsigned DIM>
void AbstractContinuumMechanicsSolver<DIM>::AllocateMatrixMemory()
{
    Vec template_vec = mrQuadMesh.GetDistributedVectorFactory()->CreateVec(mProblemDimension);

    ///////////////////////////
    // three vectors
    ///////////////////////////
    VecDuplicate(template_vec, &mResidualVector);
    VecDuplicate(mResidualVector, &mLinearSystemRhsVector);
    // the one is only allocated if it will be needed (in ApplyDirichletBoundaryConditions),
    // depending on whether the matrix is kept symmetric.
    mDirichletBoundaryConditionsVector = nullptr;
    PetscTools::Destroy(template_vec);

    ///////////////////////////
    // two matrices
    ///////////////////////////

    int lo, hi;
    VecGetOwnershipRange(mResidualVector, &lo, &hi);
    PetscInt local_size = hi - lo;


    if (DIM==2)
    {
        // 2D: N elements around a point => 7N+3 non-zeros in that row? Assume N<=10 (structured mesh would have N_max=6) => 73.
        unsigned num_non_zeros = std::min(75u, mNumDofs);

        PetscTools::SetupMat(mSystemLhsMatrix, mNumDofs, mNumDofs, num_non_zeros, local_size, local_size);
        PetscTools::SetupMat(mPreconditionMatrix, mNumDofs, mNumDofs, num_non_zeros, local_size, local_size);
    }
    else
    {
        assert(DIM==3);

        // in 3d we get the number of containing elements for each node and use that to obtain an upper bound
        // for the number of non-zeros for each DOF associated with that node.

        int* num_non_zeros_each_row = new int[mNumDofs];
        for (unsigned i=0; i<mNumDofs; i++)
        {
            num_non_zeros_each_row[i] = 0;
        }

        for (typename AbstractMesh<DIM,DIM>::NodeIterator iter = mrQuadMesh.GetNodeIteratorBegin();
             iter != mrQuadMesh.GetNodeIteratorEnd();
             ++iter)
        {
            // this upper bound neglects the fact that two containing elements will share the same nodes..
            // 4 = max num dofs associated with this node
            // 30 = 3*9+3 = 3 dimensions x 9 other nodes on this element   +  3 vertices with a pressure unknown
            unsigned num_non_zeros_upper_bound = 4 + 30*iter->GetNumContainingElements();

            num_non_zeros_upper_bound = std::min(num_non_zeros_upper_bound, mNumDofs);

            unsigned i = iter->GetIndex();

            num_non_zeros_each_row[mProblemDimension*i + 0] = num_non_zeros_upper_bound;
            num_non_zeros_each_row[mProblemDimension*i + 1] = num_non_zeros_upper_bound;
            num_non_zeros_each_row[mProblemDimension*i + 2] = num_non_zeros_upper_bound;

            if (mCompressibilityType==INCOMPRESSIBLE)
            {
                if (!iter->IsInternal())
                {
                    num_non_zeros_each_row[mProblemDimension*i + 3] = num_non_zeros_upper_bound;
                }
                else
                {
                    num_non_zeros_each_row[mProblemDimension*i + 3] = 1;
                }
            }
        }

        // NOTE: PetscTools::SetupMat() or the below creates a MATAIJ matrix, which means the matrix will
        // be of type MATSEQAIJ if num_procs=1 and MATMPIAIJ otherwise. In the former case
        // MatSeqAIJSetPreallocation MUST be called [MatMPIAIJSetPreallocation will have
        // no effect (silently)], and vice versa in the latter case

        /// We want to allocate different numbers of non-zeros per row, which means
        /// PetscTools::SetupMat isn't that useful. We could call
        //PetscTools::SetupMat(mSystemLhsMatrix, mNumDofs, mNumDofs, 0, PETSC_DECIDE, PETSC_DECIDE);
        //PetscTools::SetupMat(mPreconditionMatrix, mNumDofs, mNumDofs, 0, PETSC_DECIDE, PETSC_DECIDE);
        /// but we would get warnings due to the lack allocation

        // possible todo: create a PetscTools::SetupMatNoAllocation()


#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
        MatCreate(PETSC_COMM_WORLD,local_size,local_size,mNumDofs,mNumDofs,&mSystemLhsMatrix);
        MatCreate(PETSC_COMM_WORLD,local_size,local_size,mNumDofs,mNumDofs,&mPreconditionMatrix);
#else //New API
        MatCreate(PETSC_COMM_WORLD,&mSystemLhsMatrix);
        MatCreate(PETSC_COMM_WORLD,&mPreconditionMatrix);
        MatSetSizes(mSystemLhsMatrix,local_size,local_size,mNumDofs,mNumDofs);
        MatSetSizes(mPreconditionMatrix,local_size,local_size,mNumDofs,mNumDofs);
#endif

        if (PetscTools::IsSequential())
        {
            MatSetType(mSystemLhsMatrix, MATSEQAIJ);
            MatSetType(mPreconditionMatrix, MATSEQAIJ);
            MatSeqAIJSetPreallocation(mSystemLhsMatrix,    PETSC_DEFAULT, num_non_zeros_each_row);
            MatSeqAIJSetPreallocation(mPreconditionMatrix, PETSC_DEFAULT, num_non_zeros_each_row);
        }
        else
        {
            int* num_non_zeros_each_row_in_diag = new int[local_size];
            int* num_non_zeros_each_row_off_diag = new int[local_size];
            for (unsigned i=0; i<unsigned(local_size); i++)
            {
                num_non_zeros_each_row_in_diag[i] = num_non_zeros_each_row[lo+i];
                num_non_zeros_each_row_off_diag[i] = num_non_zeros_each_row[lo+i];
                // In the on process ("diagonal block") there cannot be more non-zero columns specified than there are rows
                if (num_non_zeros_each_row_in_diag[i] > local_size)
                {
                    num_non_zeros_each_row_in_diag[i] = local_size;
                }
            }

            MatSetType(mSystemLhsMatrix, MATMPIAIJ);
            MatSetType(mPreconditionMatrix, MATMPIAIJ);
            MatMPIAIJSetPreallocation(mSystemLhsMatrix,    PETSC_DEFAULT, num_non_zeros_each_row_in_diag, PETSC_DEFAULT, num_non_zeros_each_row_off_diag);
            MatMPIAIJSetPreallocation(mPreconditionMatrix, PETSC_DEFAULT, num_non_zeros_each_row_in_diag, PETSC_DEFAULT, num_non_zeros_each_row_off_diag);
        }

        MatSetFromOptions(mSystemLhsMatrix);
        MatSetFromOptions(mPreconditionMatrix);
#if (PETSC_VERSION_MAJOR == 3) //PETSc 3.x.x
        MatSetOption(mSystemLhsMatrix, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE);
        MatSetOption(mPreconditionMatrix, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE);
#else
        MatSetOption(mSystemLhsMatrix, MAT_IGNORE_OFF_PROC_ENTRIES);
        MatSetOption(mPreconditionMatrix, MAT_IGNORE_OFF_PROC_ENTRIES);
#endif

        //unsigned total_non_zeros = 0;
        //for (unsigned i=0; i<mNumDofs; i++)
        //{
        //   total_non_zeros += num_non_zeros_each_row[i];
        //}
        //std::cout << total_non_zeros << " versus " << 500*mNumDofs << "\n" << std::flush;

        delete [] num_non_zeros_each_row;
    }
}
#endif // ABSTRACTCONTINUUMMECHANICSSOLVER_HPP_
