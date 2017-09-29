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

#ifndef ELLIPTICBOXDOMAINPDESYSTEMMODIFIER_HPP_
#define ELLIPTICBOXDOMAINPDESYSTEMMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractBoxDomainPdeSystemModifier.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "PetscTools.hpp"
#include "FileFinder.hpp"
#include "LinearEllipticPdeSystemSolver.hpp"

#include "Debug.hpp"

/**
 * A modifier class in which a linear elliptic PDE coupled to a cell-based simulation
 * is solved on a coarse domain.
 *
 * The finite element mesh used to solve the PDE numerically is a fixed tessellation of
 * a cuboid (box), which must be supplied to the constructor. The value of the dependent
 * variable is interpolated between coarse mesh nodes to obtain a value at each cell,
 * which is stored and updated in a CellData item.
 *
 * At each time step the boundary condition supplied to the constructor may be imposed
 * either on the boundary of the box domain, or on the boundary of the cell population
 * (which is assumed to lie within the box domain). This choice can be made using the
 * AbstractBoxDomainPdeSystemModifier method SetBcsOnBoxBoundary(), which is inherited by this
 * class.
 *
 * Examples of PDEs in the source folder that can be solved using this class are
 * AveragedSourceEllipticPde, VolumeDependentAveragedSourceEllipticPde and UniformSourceEllipticPde.
 */
template<unsigned DIM, unsigned PROBLEM_DIM=1>
class EllipticBoxDomainPdeSystemModifier : public AbstractBoxDomainPdeSystemModifier<DIM, PROBLEM_DIM>
{
    friend class TestEllipticBoxDomainPdeSystemModifier;

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
        archive & boost::serialization::base_object<AbstractBoxDomainPdeSystemModifier<DIM, PROBLEM_DIM> >(*this);
MARK;
    }

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
    EllipticBoxDomainPdeSystemModifier(
        boost::shared_ptr<AbstractLinearPdeSystem<DIM,DIM,PROBLEM_DIM> > pPdeSystem=boost::shared_ptr<AbstractLinearPdeSystem<DIM,DIM,PROBLEM_DIM> >(),
        std::vector<boost::shared_ptr<AbstractBoundaryCondition<DIM> > > pBoundaryConditions=std::vector<boost::shared_ptr<AbstractBoundaryCondition<DIM> > >(),
        bool isNeumannBoundaryCondition=true,
        boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid=boost::shared_ptr<ChasteCuboid<DIM> >(),
        double stepSize=1.0,
        Vec solution=nullptr);

    /**
     * Destructor.
     */
    virtual ~EllipticBoxDomainPdeSystemModifier();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specifies what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specifies what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Helper method to construct the boundary conditions container for the PDE.
     *
     * @param rCellPopulation reference to the cell population
     *
     * @return the full boundary conditions container
     */
    virtual std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,PROBLEM_DIM> > ConstructBoundaryConditionsContainer(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

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
EllipticBoxDomainPdeSystemModifier<DIM, PROBLEM_DIM>::EllipticBoxDomainPdeSystemModifier(
    boost::shared_ptr<AbstractLinearPdeSystem<DIM, DIM, PROBLEM_DIM> > pPdeSystem,
    std::vector<boost::shared_ptr<AbstractBoundaryCondition<DIM> > > pBoundaryConditions,
    bool isNeumannBoundaryCondition,
    boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid,
    double stepSize,
    Vec solution)
    : AbstractBoxDomainPdeSystemModifier<DIM, PROBLEM_DIM>(pPdeSystem,
                                                           pBoundaryConditions,
                                                           isNeumannBoundaryCondition,
                                                           pMeshCuboid,
                                                           stepSize,
                                                           solution)
{
MARK;
}

template <unsigned DIM, unsigned PROBLEM_DIM>
EllipticBoxDomainPdeSystemModifier<DIM,PROBLEM_DIM>::~EllipticBoxDomainPdeSystemModifier()
{
}

template<unsigned DIM, unsigned PROBLEM_DIM>
void EllipticBoxDomainPdeSystemModifier<DIM, PROBLEM_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    // Set up boundary conditions
    std::shared_ptr<BoundaryConditionsContainer<DIM, DIM, PROBLEM_DIM> > p_bcc = ConstructBoundaryConditionsContainer(rCellPopulation);

    this->UpdateCellPdeElementMap(rCellPopulation);

    // When using a PDE mesh which doesn't coincide with the cells, we must set up the source terms before solving the PDE.
    // Pass in already updated CellPdeElementMap to speed up finding cells.
    this->SetUpSourceTermsForAveragedSourcePde(this->mpFeMesh, &this->mCellPdeElementMap);

    // Use LinearEllipticPdeSystemSolver as Averaged Source PDE
    ///\todo allow other PDE classes to be used with this modifier
    LinearEllipticPdeSystemSolver<DIM, DIM, PROBLEM_DIM> solver(this->mpFeMesh,
                                                                boost::static_pointer_cast<AbstractLinearEllipticPdeSystem<DIM,DIM,PROBLEM_DIM> >(this->GetPdeSystem()).get(),
                                                                p_bcc.get());

    ///\todo Use solution at previous time step as an initial guess for Solve()
    Vec old_solution_copy = this->mSolution;
    this->mSolution = solver.Solve();
    if (old_solution_copy != nullptr)
    {
        PetscTools::Destroy(old_solution_copy);
    }

    this->UpdateCellData(rCellPopulation);
}

template<unsigned DIM, unsigned PROBLEM_DIM>
void EllipticBoxDomainPdeSystemModifier<DIM, PROBLEM_DIM>::SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory)
{
    AbstractBoxDomainPdeSystemModifier<DIM, PROBLEM_DIM>::SetupSolve(rCellPopulation,outputDirectory);

    // Call these  methods to solve the PDE on the initial step and output the results
    UpdateAtEndOfTimeStep(rCellPopulation);
    this->UpdateAtEndOfOutputTimeStep(rCellPopulation);
}

template<unsigned DIM, unsigned PROBLEM_DIM>
std::shared_ptr<BoundaryConditionsContainer<DIM, DIM, PROBLEM_DIM> > EllipticBoxDomainPdeSystemModifier<DIM, PROBLEM_DIM>::ConstructBoundaryConditionsContainer(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    std::shared_ptr<BoundaryConditionsContainer<DIM, DIM, PROBLEM_DIM> > p_bcc(new BoundaryConditionsContainer<DIM,DIM,PROBLEM_DIM>(false));

    // To be well-defined, elliptic PDE problems on box domains require at least some Dirichlet boundary conditions
    ///\todo Replace this assertion with an exception in the constructor
    assert(!(this->IsNeumannBoundaryCondition()));

    if (!this->mSetBcsOnBoxBoundary)
    {
        // Get the set of coarse element indices that contain cells
        std::set<unsigned> coarse_element_indices_in_map;
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            coarse_element_indices_in_map.insert(this->mCellPdeElementMap[*cell_iter]);
        }

        // Find the node indices associated with elements whose indices are NOT in the set coarse_element_indices_in_map
        std::set<unsigned> coarse_mesh_boundary_node_indices;
        for (unsigned i=0; i<this->mpFeMesh->GetNumElements(); i++)
        {
            if (coarse_element_indices_in_map.find(i) == coarse_element_indices_in_map.end())
            {
                Element<DIM,DIM>* p_element = this->mpFeMesh->GetElement(i);
                for (unsigned j=0; j<DIM+1; j++)
                {
                    unsigned node_index = p_element->GetNodeGlobalIndex(j);
                    coarse_mesh_boundary_node_indices.insert(node_index);
                }
            }
        }

        // Apply boundary condition to the nodes in the set coarse_mesh_boundary_node_indices
        for (std::set<unsigned>::iterator iter = coarse_mesh_boundary_node_indices.begin();
             iter != coarse_mesh_boundary_node_indices.end();
             ++iter)
        {
            // Loop over PDEs
            for (unsigned pde_index=0; pde_index<PROBLEM_DIM; pde_index++)
            {
                p_bcc->AddDirichletBoundaryCondition(this->mpFeMesh->GetNode(*iter), this->mpBoundaryConditions[pde_index].get(), pde_index, false);
            }
        }
    }
    else // Apply BC at boundary nodes of box domain FE mesh
    {
        for (typename TetrahedralMesh<DIM, DIM>::BoundaryNodeIterator node_iter = this->mpFeMesh->GetBoundaryNodeIteratorBegin();
             node_iter != this->mpFeMesh->GetBoundaryNodeIteratorEnd();
             ++node_iter)
        {
            // Loop over PDEs
            for (unsigned pde_index=0; pde_index<PROBLEM_DIM; pde_index++)
            {
                p_bcc->AddDirichletBoundaryCondition(*node_iter, this->mpBoundaryConditions[pde_index].get(), pde_index);
            }
        }
    }

    return p_bcc;
}

template<unsigned DIM, unsigned PROBLEM_DIM>
void EllipticBoxDomainPdeSystemModifier<DIM, PROBLEM_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractBoxDomainPdeSystemModifier<DIM, PROBLEM_DIM>::OutputSimulationModifierParameters(rParamsFile);
}

#endif /*ELLIPTICBOXDOMAINPDESYSTEMMODIFIER_HPP_*/
