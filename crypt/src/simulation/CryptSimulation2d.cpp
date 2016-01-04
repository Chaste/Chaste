/*

Copyright (c) 2005-2016, University of Oxford.
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

#include "CryptSimulation2d.hpp"
#include "WntConcentration.hpp"
#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisOne.hpp"
#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo.hpp"
#include "CellBetaCateninWriter.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "SmartPointers.hpp"

CryptSimulation2d::CryptSimulation2d(AbstractCellPopulation<2>& rCellPopulation,
                                     bool deleteCellPopulationInDestructor,
                                     bool initialiseCells)
    : OffLatticeSimulation<2>(rCellPopulation,
                             deleteCellPopulationInDestructor,
                             initialiseCells),
      mUsingMeshBasedCellPopulation(false)
{
    /* Throw an exception message if not using a  MeshBasedCellPopulation or a VertexBasedCellPopulation.
     * This is to catch NodeBasedCellPopulations as AbstactOnLatticeBasedCellPopulations are caught in
     * the OffLatticeSimulation constructor.
     */
    if ( (dynamic_cast<VertexBasedCellPopulation<2>*>(&rCellPopulation) == NULL)
         &&(dynamic_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation) == NULL) )
    {
        EXCEPTION("CryptSimulation2d is to be used with MeshBasedCellPopulation or VertexBasedCellPopulation (or subclasses) only");
    }

    if (dynamic_cast<MeshBasedCellPopulation<2>*>(&mrCellPopulation))
    {
        mUsingMeshBasedCellPopulation = true;
    }

    if (!mDeleteCellPopulationInDestructor)
    {
        // Pass a CryptSimulationBoundaryCondition object into mBoundaryConditions
        MAKE_PTR_ARGS(CryptSimulationBoundaryCondition<2>, p_bc, (&rCellPopulation));
        AddCellPopulationBoundaryCondition(p_bc);
    }
}

CryptSimulation2d::~CryptSimulation2d()
{
}

c_vector<double, 2> CryptSimulation2d::CalculateCellDivisionVector(CellPtr pParentCell)
{
    if (mUsingMeshBasedCellPopulation)
    {
        // Location of parent and daughter cells
        c_vector<double, 2> parent_coords = mrCellPopulation.GetLocationOfCellCentre(pParentCell);
        c_vector<double, 2> daughter_coords;

        // Get separation parameter
        double separation =
            static_cast<MeshBasedCellPopulation<2>*>(&mrCellPopulation)->GetMeinekeDivisionSeparation();

        // Make a random direction vector of the required length
        c_vector<double, 2> random_vector;

        /*
         * Pick a random direction and move the parent cell backwards by 0.5*separation
         * in that direction and return the position of the daughter cell 0.5*separation
         * forwards in that direction.
         */
        double random_angle = RandomNumberGenerator::Instance()->ranf();
        random_angle *= 2.0*M_PI;

        random_vector(0) = 0.5*separation*cos(random_angle);
        random_vector(1) = 0.5*separation*sin(random_angle);

        c_vector<double, 2> proposed_new_parent_coords = parent_coords - random_vector;
        c_vector<double, 2> proposed_new_daughter_coords = parent_coords + random_vector;

        if ((proposed_new_parent_coords(1) >= 0.0) && (proposed_new_daughter_coords(1) >= 0.0))
        {
            // We are not too close to the bottom of the cell population, so move parent
            parent_coords = proposed_new_parent_coords;
            daughter_coords = proposed_new_daughter_coords;
        }
        else
        {
            proposed_new_daughter_coords = parent_coords + 2.0*random_vector;
            while (proposed_new_daughter_coords(1) < 0.0)
            {
                random_angle = RandomNumberGenerator::Instance()->ranf();
                random_angle *= 2.0*M_PI;

                random_vector(0) = separation*cos(random_angle);
                random_vector(1) = separation*sin(random_angle);
                proposed_new_daughter_coords = parent_coords + random_vector;
            }
            daughter_coords = proposed_new_daughter_coords;
        }

        assert(daughter_coords(1) >= 0.0); // to make sure dividing cells stay in the cell population
        assert(parent_coords(1) >= 0.0);   // to make sure dividing cells stay in the cell population

        // Set the parent to use this location
        ChastePoint<2> parent_coords_point(parent_coords);

        unsigned node_index = mrCellPopulation.GetLocationIndexUsingCell(pParentCell);
        mrCellPopulation.SetNode(node_index, parent_coords_point);

        return daughter_coords;
    }
    else // using a VertexBasedCellPopulation
    {
        // Let's check we're a VertexBasedCellPopulation when we're in debug mode...
        assert(dynamic_cast<VertexBasedCellPopulation<2>*>(&(this->mrCellPopulation)));

        VertexBasedCellPopulation<2>* p_vertex_population = dynamic_cast<VertexBasedCellPopulation<2>*>(&(this->mrCellPopulation));
        c_vector<double, 2> axis_of_division = p_vertex_population->
                GetVertexBasedDivisionRule()->CalculateCellDivisionVector(pParentCell, *p_vertex_population);

        // We don't need to prescribe how 'stem' cells divide if Wnt is present
        bool is_wnt_included = WntConcentration<2>::Instance()->IsWntSetUp();
        if (!is_wnt_included)
        {
            WntConcentration<2>::Destroy();
            if (pParentCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
            {
                axis_of_division(0) = 1.0;
                axis_of_division(1) = 0.0;
            }
        }
        return axis_of_division;
    }
}

void CryptSimulation2d::SetupSolve()
{
    // First call method on base class
    OffLatticeSimulation<2>::SetupSolve();

    /*
     * To check if beta-catenin results will be written to file, we test if the first
     * cell has a cell-cycle model that is a subclass of AbstractVanLeeuwen2009WntSwatCellCycleModel.
     * In doing so, we assume that all cells in the simulation have the same cell-cycle
     * model.
     */
    if (dynamic_cast<AbstractVanLeeuwen2009WntSwatCellCycleModel*>(this->mrCellPopulation.Begin()->GetCellCycleModel()))
    {
        mrCellPopulation.AddCellWriter<CellBetaCateninWriter>();
    }

    if (dynamic_cast<AbstractVanLeeuwen2009WntSwatCellCycleModel*>(mrCellPopulation.Begin()->GetCellCycleModel()))
    {
        *mpVizSetupFile << "BetaCatenin\n";
    }
}

void CryptSimulation2d::UseJiggledBottomCells()
{
    // The CryptSimulationBoundaryCondition object is the first element of mBoundaryConditions
    boost::static_pointer_cast<CryptSimulationBoundaryCondition<2> >(mBoundaryConditions[0])->SetUseJiggledBottomCells(true);
}

void CryptSimulation2d::SetBottomCellAncestors()
{
    /*
     * We use a different height threshold depending on which type of cell
     * population we are using, a MeshBasedCellPopulationWithGhostNodes or
     * a VertexBasedCellPopulation.
     */
    double threshold_height = 1.0;
    if (mUsingMeshBasedCellPopulation)
    {
        threshold_height = 0.5;
    }

    unsigned index = 0;
    for (AbstractCellPopulation<2>::Iterator cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {
        if (mrCellPopulation.GetLocationOfCellCentre(*cell_iter)[1] < threshold_height)
        {
            MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (index++));
            cell_iter->SetAncestor(p_cell_ancestor);
        }
    }
}

void CryptSimulation2d::OutputSimulationParameters(out_stream& rParamsFile)
{
    double width = mrCellPopulation.GetWidth(0);
    bool use_jiggled_bottom_cells = boost::static_pointer_cast<CryptSimulationBoundaryCondition<2> >(mBoundaryConditions[0])->GetUseJiggledBottomCells();

    *rParamsFile << "\t\t<CryptCircumference>" << width << "</CryptCircumference>\n";
    *rParamsFile << "\t\t<UseJiggledBottomCells>" << use_jiggled_bottom_cells << "</UseJiggledBottomCells>\n";

    // Call method on direct parent class
    OffLatticeSimulation<2>::OutputSimulationParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CryptSimulation2d)
