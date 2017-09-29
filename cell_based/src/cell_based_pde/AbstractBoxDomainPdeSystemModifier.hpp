/*

Copyright (c) 2005-2017, University of Oxford.
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

#ifndef ABSTRACTBOXDOMAINPDESYSTEMMODIFIER_HPP_
#define ABSTRACTBOXDOMAINPDESYSTEMMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "AbstractPdeSystemModifier.hpp"
#include "ReplicatableVector.hpp"
#include "LinearBasisFunction.hpp"

/**
 * An abstract modifier class containing functionality common to EllipticBoxDomainPdeSystemModifier
 * and ParabolicBoxDomainPdeSystemModifier, which both solve a linear elliptic or parabolic PDE
 * coupled to a cell-based simulation on a coarse domain.
 */
template<unsigned DIM, unsigned PROBLEM_DIM=1>
class AbstractBoxDomainPdeSystemModifier : public AbstractPdeSystemModifier<DIM, PROBLEM_DIM>
{
    friend class TestEllipticBoxDomainPdeSystemModifier;
    friend class TestParabolicBoxDomainPdeSystemModifier;
    friend class TestOffLatticeSimulationWithPdes;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
MARK;
        archive & boost::serialization::base_object<AbstractPdeSystemModifier<DIM, PROBLEM_DIM> >(*this);
MARK;
//        archive & mCellPdeElementMap; ///\todo #2930
        archive & mpMeshCuboid;
        archive & mStepSize;
        archive & mSetBcsOnBoxBoundary;
MARK;
    }

protected:

    /** Map between cells and the elements of the FE mesh containing them. */
    std::map<CellPtr, unsigned> mCellPdeElementMap;

    /**
     * Pointer to a ChasteCuboid storing the outer boundary for the FE mesh.
     * Stored as a member to facilitate archiving.
     */
    boost::shared_ptr<ChasteCuboid<DIM> > mpMeshCuboid;

    /**
     * The step size to be used in the FE mesh.
     * Stored as a member to facilitate archiving.
     */
    double mStepSize;

    /**
     * Whether to set the boundary condition on the edge of the box domain rather than the cell population.
     * Default to true.
     */
    bool mSetBcsOnBoxBoundary;

public:

    /**
     * Constructor.
     *
     * @param pPdeSystem A shared pointer to a linear PDE system object (defaults to NULL)
     * @param pBoundaryConditions A vector of shared pointers to abstract boundary conditions
     *     (defaults to NULL, corresponding to a constant boundary condition with value zero)
     * @param isNeumannBoundaryCondition Whether the boundary condition is Neumann (defaults to true)
     * @param pMeshCuboid A shared pointer to a ChasteCuboid specifying the outer boundary for the FE mesh (defaults to NULL)
     * @param stepSize step size to be used in the FE mesh (defaults to 1.0, i.e. the default cell size)
     * @param solution solution vector (defaults to NULL)
     */
    AbstractBoxDomainPdeSystemModifier(
        boost::shared_ptr<AbstractLinearPdeSystem<DIM,DIM,PROBLEM_DIM> > pPdeSystem=boost::shared_ptr<AbstractLinearPdeSystem<DIM,DIM,PROBLEM_DIM> >(),
        std::vector<boost::shared_ptr<AbstractBoundaryCondition<DIM> > > pBoundaryConditions=std::vector<boost::shared_ptr<AbstractBoundaryCondition<DIM> > >(),
        bool isNeumannBoundaryCondition=true,
        boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid=boost::shared_ptr<ChasteCuboid<DIM> >(),
        double stepSize=1.0,
        Vec solution=nullptr);

    /**
     * Destructor.
     */
    virtual ~AbstractBoxDomainPdeSystemModifier();

    /**
     * @return mStepSize.
     */
    double GetStepSize();

    /**
     * Set mSetBcsOnCoarseBoundary.
     *
     * @param setBcsOnBoxBoundary whether to set the boundary condition on the edge of the box domain rather than the cell population
     */
    void SetBcsOnBoxBoundary(bool setBcsOnBoxBoundary);

    /**
     * @return mSetBcsOnCoarseBoundary.
     */
    bool AreBcsSetOnBoxBoundary();

    /**
     * Overridden SetupSolve() method.
     *
     * Specifies what to do in the simulation before the start of the time loop.
     *
     * Here we just initialize the Cell PDE element map
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Helper method to generate the mesh.
     *
     * @param pMeshCuboid the outer boundary for the FE mesh.
     * @param stepSize the step size to be used in the FE mesh.
     */
    void GenerateFeMesh(boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid, double stepSize);

    /**
     * Helper method to copy the PDE solution to CellData
     *
     * Here we need to interpolate from the FE mesh onto the cells.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Initialise mCellPdeElementMap.
     *
     * @param rCellPopulation reference to the cell population
     */
    void InitialiseCellPdeElementMap(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Update the mCellPdeElementMap
     *
     * This method should be called before sending the element map to a PDE class
     * to ensure map is up to date.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateCellPdeElementMap(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

/*
 * As this class is templated over PROBLEM_DIM, we put the implementation
 * in the header file to avoid explicit instantiation.
 */

template<unsigned DIM, unsigned PROBLEM_DIM>
AbstractBoxDomainPdeSystemModifier<DIM, PROBLEM_DIM>::AbstractBoxDomainPdeSystemModifier(
    boost::shared_ptr<AbstractLinearPdeSystem<DIM, DIM, PROBLEM_DIM> > pPdeSystem,
    std::vector<boost::shared_ptr<AbstractBoundaryCondition<DIM> > > pBoundaryConditions,
    bool isNeumannBoundaryCondition,
    boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid,
    double stepSize,
    Vec solution)
    : AbstractPdeSystemModifier<DIM, PROBLEM_DIM>(pPdeSystem,
                                                  pBoundaryConditions,
                                                  isNeumannBoundaryCondition,
                                                  solution),
      mpMeshCuboid(pMeshCuboid),
      mStepSize(stepSize),
      mSetBcsOnBoxBoundary(true)
{
MARK;
    if (pMeshCuboid)
    {
        // We only need to generate mpFeMesh once, as it does not vary with time
        this->GenerateFeMesh(mpMeshCuboid, mStepSize);
        this->mDeleteFeMesh = true;
    }
}

template<unsigned DIM, unsigned PROBLEM_DIM>
AbstractBoxDomainPdeSystemModifier<DIM, PROBLEM_DIM>::~AbstractBoxDomainPdeSystemModifier()
{
}

template<unsigned DIM, unsigned PROBLEM_DIM>
double AbstractBoxDomainPdeSystemModifier<DIM, PROBLEM_DIM>::GetStepSize()
{
     return mStepSize;
}

template<unsigned DIM, unsigned PROBLEM_DIM>
void AbstractBoxDomainPdeSystemModifier<DIM, PROBLEM_DIM>::SetBcsOnBoxBoundary(bool setBcsOnBoxBoundary)
{
    mSetBcsOnBoxBoundary = setBcsOnBoxBoundary;
}

template<unsigned DIM, unsigned PROBLEM_DIM>
bool AbstractBoxDomainPdeSystemModifier<DIM, PROBLEM_DIM>::AreBcsSetOnBoxBoundary()
{
    return mSetBcsOnBoxBoundary;
}

template<unsigned DIM, unsigned PROBLEM_DIM>
void AbstractBoxDomainPdeSystemModifier<DIM, PROBLEM_DIM>::SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory)
{
    AbstractPdeSystemModifier<DIM,PROBLEM_DIM>::SetupSolve(rCellPopulation, outputDirectory);
    InitialiseCellPdeElementMap(rCellPopulation);
}

template<unsigned DIM, unsigned PROBLEM_DIM>
void AbstractBoxDomainPdeSystemModifier<DIM, PROBLEM_DIM>::GenerateFeMesh(boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid, double stepSize)
{
    // Create a regular coarse tetrahedral mesh
    this->mpFeMesh = new TetrahedralMesh<DIM,DIM>();
    switch (DIM)
    {
        case 1:
            this->mpFeMesh->ConstructRegularSlabMesh(stepSize, pMeshCuboid->GetWidth(0));
            break;
        case 2:
            this->mpFeMesh->ConstructRegularSlabMesh(stepSize, pMeshCuboid->GetWidth(0), pMeshCuboid->GetWidth(1));
            break;
        case 3:
            this->mpFeMesh->ConstructRegularSlabMesh(stepSize, pMeshCuboid->GetWidth(0), pMeshCuboid->GetWidth(1), pMeshCuboid->GetWidth(2));
            break;
        default:
            NEVER_REACHED;
    }

    // Get centroid of meshCuboid
    ChastePoint<DIM> upper = pMeshCuboid->rGetUpperCorner();
    ChastePoint<DIM> lower = pMeshCuboid->rGetLowerCorner();
    c_vector<double,DIM> centre_of_cuboid = 0.5*(upper.rGetLocation() + lower.rGetLocation());

    // Find the centre of the PDE mesh
    c_vector<double,DIM> centre_of_coarse_mesh = zero_vector<double>(DIM);
    for (unsigned i=0; i<this->mpFeMesh->GetNumNodes(); i++)
    {
        centre_of_coarse_mesh += this->mpFeMesh->GetNode(i)->rGetLocation();
    }
    centre_of_coarse_mesh /= this->mpFeMesh->GetNumNodes();

    // Now move the mesh to the correct location
    this->mpFeMesh->Translate(centre_of_cuboid - centre_of_coarse_mesh);
}

template<unsigned DIM, unsigned PROBLEM_DIM>
void AbstractBoxDomainPdeSystemModifier<DIM, PROBLEM_DIM>::UpdateCellData(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    // Store the PDE solution in an accessible form
    ReplicatableVector solution_repl(this->mSolution);

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // The cells are not nodes of the mesh, so we must interpolate
        double solution_at_cell = 0.0;

        // Find the element in the FE mesh that contains this cell. CellElementMap has been updated so use this.
        unsigned elem_index = mCellPdeElementMap[*cell_iter];
        Element<DIM,DIM>* p_element = this->mpFeMesh->GetElement(elem_index);

        const ChastePoint<DIM>& node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

        c_vector<double,DIM+1> weights = p_element->CalculateInterpolationWeights(node_location);

        for (unsigned i=0; i<DIM+1; i++)
        {
            double nodal_value = solution_repl[p_element->GetNodeGlobalIndex(i)];
            solution_at_cell += nodal_value * weights(i);
        }

        cell_iter->GetCellData()->SetItem(this->mDependentVariableNames[0], solution_at_cell);

        if (this->mOutputGradient)
        {
            // Now calculate the gradient of the solution and store this in CellVecData
            c_vector<double, DIM> solution_gradient = zero_vector<double>(DIM);

            // Calculate the basis functions at any point (e.g. zero) in the element
            c_matrix<double, DIM, DIM> jacobian, inverse_jacobian;
            double jacobian_det;
            this->mpFeMesh->GetInverseJacobianForElement(elem_index, jacobian, jacobian_det, inverse_jacobian);
            const ChastePoint<DIM> zero_point;
            c_matrix<double, DIM, DIM+1> grad_phi;
            LinearBasisFunction<DIM>::ComputeTransformedBasisFunctionDerivatives(zero_point, inverse_jacobian, grad_phi);

            for (unsigned node_index=0; node_index<DIM+1; node_index++)
            {
                double nodal_value = solution_repl[p_element->GetNodeGlobalIndex(node_index)];

                for (unsigned j=0; j<DIM; j++)
                {
                    solution_gradient(j) += nodal_value* grad_phi(j, node_index);
                }
            }

            switch (DIM)
            {
                case 1:
                    cell_iter->GetCellData()->SetItem(this->mDependentVariableNames[0]+"_grad_x", solution_gradient(0));
                    break;
                case 2:
                    cell_iter->GetCellData()->SetItem(this->mDependentVariableNames[0]+"_grad_x", solution_gradient(0));
                    cell_iter->GetCellData()->SetItem(this->mDependentVariableNames[0]+"_grad_y", solution_gradient(1));
                    break;
                case 3:
                    cell_iter->GetCellData()->SetItem(this->mDependentVariableNames[0]+"_grad_x", solution_gradient(0));
                    cell_iter->GetCellData()->SetItem(this->mDependentVariableNames[0]+"_grad_y", solution_gradient(1));
                    cell_iter->GetCellData()->SetItem(this->mDependentVariableNames[0]+"_grad_z", solution_gradient(2));
                    break;
                default:
                    NEVER_REACHED;
            }
        }
    }
}

template<unsigned DIM, unsigned PROBLEM_DIM>
void AbstractBoxDomainPdeSystemModifier<DIM, PROBLEM_DIM>::InitialiseCellPdeElementMap(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    mCellPdeElementMap.clear();

    // Find the element of mpFeMesh that contains each cell and populate mCellPdeElementMap
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        const ChastePoint<DIM>& r_position_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
        unsigned elem_index = this->mpFeMesh->GetContainingElementIndex(r_position_of_cell);
        mCellPdeElementMap[*cell_iter] = elem_index;
    }
}

template<unsigned DIM, unsigned PROBLEM_DIM>
void AbstractBoxDomainPdeSystemModifier<DIM, PROBLEM_DIM>::UpdateCellPdeElementMap(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    // Find the element of mpCoarsePdeMesh that contains each cell and populate mCellPdeElementMap
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        const ChastePoint<DIM>& r_position_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
        unsigned elem_index = this->mpFeMesh->GetContainingElementIndexWithInitialGuess(r_position_of_cell, mCellPdeElementMap[*cell_iter]);
        mCellPdeElementMap[*cell_iter] = elem_index;
    }
}

template<unsigned DIM, unsigned PROBLEM_DIM>
void AbstractBoxDomainPdeSystemModifier<DIM, PROBLEM_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractPdeSystemModifier<DIM,PROBLEM_DIM>::OutputSimulationModifierParameters(rParamsFile);
}

#endif /*ABSTRACTBOXDOMAINPDESYSTEMMODIFIER_HPP_*/
