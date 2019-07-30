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

#include "CryptSimulation2d.hpp"
#include "CellAncestor.hpp"
#include "CellBetaCateninWriter.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "SmartPointers.hpp"
#include "StemCellProliferativeType.hpp"
#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisOne.hpp"
#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo.hpp"
#include "WntConcentration.hpp"

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
    if ((dynamic_cast<VertexBasedCellPopulation<2>*>(&rCellPopulation) == nullptr)
        && (dynamic_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation) == nullptr))
    {
        EXCEPTION("CryptSimulation2d is to be used with MeshBasedCellPopulation or VertexBasedCellPopulation (or subclasses) only");
    }

    if (dynamic_cast<MeshBasedCellPopulation<2>*>(&mrCellPopulation))
    {
        mUsingMeshBasedCellPopulation = true;

        MAKE_PTR(CryptCentreBasedDivisionRule<2>, p_centre_div_rule);
        static_cast<MeshBasedCellPopulation<2>*>(&mrCellPopulation)->SetCentreBasedDivisionRule(p_centre_div_rule);
    }
    else // VertexBasedCellPopulation
    {
        MAKE_PTR(CryptVertexBasedDivisionRule<2>, p_vertex_div_rule);
        static_cast<VertexBasedCellPopulation<2>*>(&mrCellPopulation)->SetVertexBasedDivisionRule(p_vertex_div_rule);
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
     *
     * \todo Make this threshold height a member variable and set it in the constructor,
     *       depending on the cell population type; this would allow us to remove mUsingMeshBasedCellPopulation
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
