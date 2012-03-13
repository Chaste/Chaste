/*

Copyright (c) 2005-2012, University of Oxford.
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

#include "OffLatticeSimulation.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

#include "Cylindrical2dMesh.hpp"
#include "Cylindrical2dVertexMesh.hpp"

#include "AbstractTwoBodyInteractionForce.hpp"
#include "CellBasedEventHandler.hpp"
#include "LogFile.hpp"
#include "Version.hpp"
#include "ExecutableSupport.hpp"

template<unsigned DIM>
OffLatticeSimulation<DIM>::OffLatticeSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
                                                bool deleteCellPopulationInDestructor,
                                                bool initialiseCells)
    : AbstractCellBasedSimulation<DIM>(rCellPopulation, deleteCellPopulationInDestructor, initialiseCells),
      mOutputNodeVelocities(false)
{
    if (!dynamic_cast<AbstractOffLatticeCellPopulation<DIM>*>(&rCellPopulation))
    {
        EXCEPTION("OffLatticeSimulations require a subclass of AbstractOffLatticeCellPopulation.");
    }

    // Different time steps are used for cell-centre and vertex-based simulations
    if (dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(&rCellPopulation))
    {
        this->mDt = 1.0/120.0; // 30 seconds
    }
    else
    {
        assert (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation));
        this->mDt = 0.002; // smaller time step required for convergence/stability
    }
}

template<unsigned DIM>
void OffLatticeSimulation<DIM>::AddForce(boost::shared_ptr<AbstractForce<DIM> > pForce)
{
    mForceCollection.push_back(pForce);
}

template<unsigned DIM>
void OffLatticeSimulation<DIM>::AddCellPopulationBoundaryCondition(boost::shared_ptr<AbstractCellPopulationBoundaryCondition<DIM> > pBoundaryCondition)
{
    mBoundaryConditions.push_back(pBoundaryCondition);
}

template<unsigned DIM>
void OffLatticeSimulation<DIM>::UpdateCellLocationsAndTopology()
{
    // Calculate forces
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::FORCE);

    // Initialise a vector of forces on node
    std::vector<c_vector<double, DIM> > forces(this->mrCellPopulation.GetNumNodes(), zero_vector<double>(DIM));

/**\todo Is it faster to preallocate and have forces as a member variable? see #1890**/
//    // First set all the forces to zero
//    for (unsigned i=0; i<forces.size(); i++)
//    {
//         forces[i].clear();
//    }
//
//    // Then resize the std::vector if the number of cells has increased or decreased
//    // (note this should be done after the above zeroing)
//    unsigned num_nodes = mrCellPopulation.GetNumNodes();
//    if (num_nodes != forces.size())
//    {
//        forces.resize(num_nodes, zero_vector<double>(DIM));
//    }

    // Now add force contributions from each AbstractForce
    for (typename std::vector<boost::shared_ptr<AbstractForce<DIM> > >::iterator iter = mForceCollection.begin();
         iter != mForceCollection.end();
         ++iter)
    {
        (*iter)->AddForceContribution(forces, this->mrCellPopulation);
    }
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::FORCE);

    // Update node positions
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::POSITION);
    UpdateNodePositions(forces);
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::POSITION);
}

template<unsigned DIM>
c_vector<double, DIM> OffLatticeSimulation<DIM>::CalculateCellDivisionVector(CellPtr pParentCell)
{
    /**
     * \todo Could remove this dynamic_cast by moving the code block below into
     * AbstractCentreBasedCellPopulation::AddCell(), allowing it to be overruled by
     * this method when overridden in subclasses. See also comment on #1093.
     */
    if (dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)))
    {
        // Location of parent and daughter cells
        c_vector<double, DIM> parent_coords = this->mrCellPopulation.GetLocationOfCellCentre(pParentCell);
        c_vector<double, DIM> daughter_coords;

        // Get separation parameter
        double separation = static_cast<AbstractCentreBasedCellPopulation<DIM>*>(&(this->mrCellPopulation))->GetMeinekeDivisionSeparation();

        // Make a random direction vector of the required length
        c_vector<double, DIM> random_vector;

        /*
         * Pick a random direction and move the parent cell backwards by 0.5*separation
         * in that direction and return the position of the daughter cell 0.5*separation
         * forwards in that direction.
         */
        switch (DIM)
        {
            case 1:
            {
                double random_direction = -1.0 + 2.0*(RandomNumberGenerator::Instance()->ranf() < 0.5);

                random_vector(0) = 0.5*separation*random_direction;
                break;
            }
            case 2:
            {
                double random_angle = 2.0*M_PI*RandomNumberGenerator::Instance()->ranf();

                random_vector(0) = 0.5*separation*cos(random_angle);
                random_vector(1) = 0.5*separation*sin(random_angle);
                break;
            }
            case 3:
            {
                double random_zenith_angle = M_PI*RandomNumberGenerator::Instance()->ranf(); // phi
                double random_azimuth_angle = 2*M_PI*RandomNumberGenerator::Instance()->ranf(); // theta

                random_vector(0) = 0.5*separation*cos(random_azimuth_angle)*sin(random_zenith_angle);
                random_vector(1) = 0.5*separation*sin(random_azimuth_angle)*sin(random_zenith_angle);
                random_vector(2) = 0.5*separation*cos(random_zenith_angle);
                break;
            }
            default:
                // This can't happen
                NEVER_REACHED;
        }

        parent_coords = parent_coords - random_vector;
        daughter_coords = parent_coords + random_vector;

        // Set the parent to use this location
        ChastePoint<DIM> parent_coords_point(parent_coords);
        unsigned node_index = this->mrCellPopulation.GetLocationIndexUsingCell(pParentCell);
        this->mrCellPopulation.SetNode(node_index, parent_coords_point);

        return daughter_coords;
    }
    else
    {
        ///\todo do something for vertex models here
        return zero_vector<double>(DIM);
    }
}

template<unsigned DIM>
void OffLatticeSimulation<DIM>::WriteVisualizerSetupFile()
{
	if(PetscTools::AmMaster())
	{
		for (unsigned i=0; i<this->mForceCollection.size(); i++)
		{
			boost::shared_ptr<AbstractForce<DIM> > p_force = this->mForceCollection[i];
			if (boost::dynamic_pointer_cast<AbstractTwoBodyInteractionForce<DIM> >(p_force))
			{
				double cutoff = (boost::static_pointer_cast<AbstractTwoBodyInteractionForce<DIM> >(p_force))->GetCutOffLength();
				*(this->mpVizSetupFile) << "Cutoff\t" << cutoff << "\n";
			}
		}

		// This is a quick and dirty check to see if the mesh is periodic
		if (dynamic_cast<MeshBasedCellPopulation<DIM>*>(&this->mrCellPopulation))
		{
		   if (dynamic_cast<Cylindrical2dMesh*>(&(dynamic_cast<MeshBasedCellPopulation<DIM>*>(&(this->mrCellPopulation))->rGetMesh())))
		   {
			   *this->mpVizSetupFile << "MeshWidth\t" << this->mrCellPopulation.GetWidth(0) << "\n";
		   }
		}
		else if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&this->mrCellPopulation))
		{
		   if (dynamic_cast<Cylindrical2dVertexMesh*>(&(dynamic_cast<VertexBasedCellPopulation<DIM>*>(&(this->mrCellPopulation))->rGetMesh())))
		   {
			   *this->mpVizSetupFile << "MeshWidth\t" << this->mrCellPopulation.GetWidth(0) << "\n";
		   }
		}
	}
}

template<unsigned DIM>
void OffLatticeSimulation<DIM>::UpdateNodePositions(const std::vector< c_vector<double, DIM> >& rNodeForces)
{
    unsigned num_nodes = this->mrCellPopulation.GetNumNodes();

    /*
     * Get the previous node positions (these may be needed when applying boundary conditions,
     * e.g. in the case of immotile cells)
     */
    std::vector<c_vector<double, DIM> > old_node_locations;
    old_node_locations.reserve(num_nodes);
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        old_node_locations[node_index] = this->mrCellPopulation.GetNode(node_index)->rGetLocation();
    }

    // Update node locations
    static_cast<AbstractOffLatticeCellPopulation<DIM>*>(&(this->mrCellPopulation))->UpdateNodeLocations(rNodeForces, this->mDt);

    // Apply any boundary conditions
    for (typename std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<DIM> > >::iterator bcs_iter = mBoundaryConditions.begin();
         bcs_iter != mBoundaryConditions.end();
         ++bcs_iter)
    {
        (*bcs_iter)->ImposeBoundaryCondition(old_node_locations);
    }

    // Verify that each boundary condition is now satisfied
    for (typename std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<DIM> > >::iterator bcs_iter = mBoundaryConditions.begin();
         bcs_iter != mBoundaryConditions.end();
         ++bcs_iter)
    {
        if (!((*bcs_iter)->VerifyBoundaryCondition()))
        {
            EXCEPTION("The cell population boundary conditions are incompatible.");
        }
    }

    // Write node velocities to file if required
    if (mOutputNodeVelocities)
    {
        if (SimulationTime::Instance()->GetTimeStepsElapsed()%this->mSamplingTimestepMultiple == 0)
        {
            *mpNodeVelocitiesFile << SimulationTime::Instance()->GetTime() << "\t";
            for (unsigned node_index=0; node_index<num_nodes; node_index++)
            {
                // We should never encounter deleted nodes due to where this method is called by Solve()
                assert(!this->mrCellPopulation.GetNode(node_index)->IsDeleted());

                // Check that results should be written for this node
                bool is_real_node = true;

                if (dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(&this->mrCellPopulation))
                {
                    if (static_cast<AbstractCentreBasedCellPopulation<DIM>*>(&(this->mrCellPopulation))->IsGhostNode(node_index))
                    {
                        // If this node is a ghost node then don't record its velocity
                        is_real_node = false;
                    }
                    else
                    {
                        // We should never encounter nodes associated with dead cells due to where this method is called by Solve()
                        assert(!this->mrCellPopulation.GetCellUsingLocationIndex(node_index)->IsDead());
                    }
                }

                // Write node data to file
                if (is_real_node)
                {
                    const c_vector<double,DIM>& position = this->mrCellPopulation.GetNode(node_index)->rGetLocation();
                    double damping_constant = static_cast<AbstractOffLatticeCellPopulation<DIM>*>(&(this->mrCellPopulation))->GetDampingConstant(node_index);
                    c_vector<double, DIM> velocity = this->mDt * rNodeForces[node_index] / damping_constant;

                    *mpNodeVelocitiesFile << node_index  << " ";
                    for (unsigned i=0; i<DIM; i++)
                    {
                        *mpNodeVelocitiesFile << position[i] << " ";
                    }
                    for (unsigned i=0; i<DIM; i++)
                    {
                        *mpNodeVelocitiesFile << velocity[i] << " ";
                    }
                }
            }
            *mpNodeVelocitiesFile << "\n";
        }
    }
}

template<unsigned DIM>
void OffLatticeSimulation<DIM>::SetupSolve()
{
    if (mOutputNodeVelocities)
    {
        OutputFileHandler output_file_handler2(this->mSimulationOutputDirectory+"/", false);
        mpNodeVelocitiesFile = output_file_handler2.OpenOutputFile("nodevelocities.dat");
    }
}

template<unsigned DIM>
void OffLatticeSimulation<DIM>::UpdateAtEndOfSolve()
{
    if (mOutputNodeVelocities)
    {
        mpNodeVelocitiesFile->close();
    }
}

template<unsigned DIM>
void OffLatticeSimulation<DIM>::OutputAdditionalSimulationSetup(out_stream& rParamsFile)
{
    // Loop over forces
    *rParamsFile << "\n\t<Forces>\n";
    for (typename std::vector<boost::shared_ptr<AbstractForce<DIM> > >::iterator iter = mForceCollection.begin();
         iter != mForceCollection.end();
         ++iter)
    {
        // Output force details
        (*iter)->OutputForceInfo(rParamsFile);
    }
    *rParamsFile << "\t</Forces>\n";

    // Loop over cell population boundary conditions
    *rParamsFile << "\n\t<CellPopulationBoundaryConditions>\n";
    for (typename std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<DIM> > >::iterator iter = mBoundaryConditions.begin();
         iter != mBoundaryConditions.end();
         ++iter)
    {
        // Output cell Boundary condition details
        (*iter)->OutputCellPopulationBoundaryConditionInfo(rParamsFile);
    }
    *rParamsFile << "\t</CellPopulationBoundaryConditions>\n";
}

template<unsigned DIM>
bool OffLatticeSimulation<DIM>::GetOutputNodeVelocities()
{
    return mOutputNodeVelocities;
}

template<unsigned DIM>
void OffLatticeSimulation<DIM>::SetOutputNodeVelocities(bool outputNodeVelocities)
{
    mOutputNodeVelocities = outputNodeVelocities;
}

template<unsigned DIM>
void OffLatticeSimulation<DIM>::OutputSimulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<OutputNodeVelocities>" << mOutputNodeVelocities << "</OutputNodeVelocities>\n";

    // Call method on direct parent class
    AbstractCellBasedSimulation<DIM>::OutputSimulationParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class OffLatticeSimulation<1>;
template class OffLatticeSimulation<2>;
template class OffLatticeSimulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(OffLatticeSimulation)
