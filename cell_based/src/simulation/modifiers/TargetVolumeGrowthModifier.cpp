/* 
  
 Copyright (c) 2005-2014, University of Oxford. 
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
 	 
  
 #include "TargetVolumeGrowthModifier.hpp"
 #include "MeshBasedCellPopulation.hpp"
  
 template<unsigned DIM>
 TargetVolumeGrowthModifier<DIM>::TargetVolumeGrowthModifier() 
 	: AbstractCellBasedSimulationModifier<DIM>() 
 { 
         mMatureCellTargetVolume = 1.0;
         mMeanSlope = 10.0; // IF I TAKE OFF THE SLOPE MODIFIER
 }

 template<unsigned DIM> 
 TargetVolumeGrowthModifier<DIM>::~TargetVolumeGrowthModifier() 
 { 
 }

 template<unsigned DIM> 
 void TargetVolumeGrowthModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)  
 { 
	UpdateTargetVolumes(rCellPopulation); 
 }
  
 template<unsigned DIM> 
 void TargetVolumeGrowthModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory) 
 { 
     /* 
     * We must update CellData in SetupSolve(), otherwise it will not have been 
     * fully initialised by the time we enter the main time loop. 
      */ 
     UpdateTargetVolumes(rCellPopulation); 
 } 
  
 template<unsigned DIM> 
 void TargetVolumeGrowthModifier<DIM>::UpdateTargetVolumes(AbstractCellPopulation<DIM,DIM>& rCellPopulation) 
 	 
 { 
     // Make sure the cell population is updated 
     ///\todo: double check that this update call doesn't break anything (i.e. counting of swaps etc.) 
     rCellPopulation.Update(); 
 	 
     if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == NULL) 
     { 
     	EXCEPTION("TargetVolumeGrowthModifier is to be used with a VertexBasedCellPopulation only"); 
     } 
 	 
     VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation); 
 	 
     for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin(); 
              elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd(); 
              ++elem_iter) 
         { 
         unsigned elem_index = elem_iter->GetIndex(); 
         //double number_of_cells = rCellPopulation.GetNumRealCells();
         UpdateTargetVolumeOfCell( p_cell_population->GetCellUsingLocationIndex(elem_index));
         } 
 } 
 	 
 template<unsigned DIM> 
 void TargetVolumeGrowthModifier<DIM>::UpdateTargetVolumeOfCell(CellPtr pCell)
 { 
     // Get target volume A of a healthy cell in S, G2 or M phase 
     double cell_target_volume = mMatureCellTargetVolume; // DO NOT CHANGE

     //<3
     double cell_slope = pCell->GetCellData()->GetItem("slope");


     // if I do not want a modifier for the slope anymore (BUT can be useful for Mao2)
     //double number_of_cells = rCellPopulation.GetNumRealCells();
     //double cell_id = pCell->GetCellId();
     //double cell_slope = mMeanSlope + 0.0*(cell_id-(number_of_cells/2.0))/number_of_cells;

 	 
     double cell_age = pCell->GetAge(); 
     double g1_duration = pCell->GetCellCycleModel()->GetG1Duration(); 
     double s_duration = pCell->GetCellCycleModel()->GetSDuration();
 	 
     // If the cell is differentiated then its G1 duration is infinite 
     if (g1_duration == DBL_MAX) // don't use magic number, compare to DBL_MAX 
     { 
         // This is just for fixed cell-cycle models, need to work out how to find the g1 duration 
         g1_duration = pCell->GetCellCycleModel()->GetTransitCellG1Duration(); 
     } 
 	 
     // If apoptic cell
     if (pCell->HasCellProperty<ApoptoticCellProperty>()) 
     { 
         // Age of cell when apoptosis begins 
         if (pCell->GetStartOfApoptosisTime() - pCell->GetBirthTime() < g1_duration) 
         {
        	 cell_target_volume *= 0.5*(1 + (pCell->GetStartOfApoptosisTime() - pCell->GetBirthTime())/g1_duration); // probleme because of my first line
         
         }
  
         // The target volume of an apoptotic cell decreases linearly to zero (and past it negative) 
         cell_target_volume = cell_target_volume - 0.5*cell_target_volume/(pCell->GetApoptosisTime())*(SimulationTime::Instance()->GetTime()-pCell->GetStartOfApoptosisTime()); 
 	 
         // Don't allow a negative target volume 
         if (cell_target_volume < 0) 
         { 
             cell_target_volume = 0; 
         } 
     } 
     else 
     { 
         if (cell_age <= g1_duration*1.1) // \todo change the 1.1
         { 
             cell_target_volume = 0.5 + cell_age/cell_slope; // \todo  determine the slope
         } 
         else 
         { 
                 // At division, daughter cells inherit the cell data array from the mother cell. Here, we assign the target volume 
                 // that we want daughter cells to have to cells that we know to divide in this time step. This is a little hack 
                 // that we might want to clean up in the future. 
                 if (pCell->ReadyToDivide()) 
                 { 
                	 cell_target_volume = 0.5*mMatureCellTargetVolume;
                 } 

         } 
         if (pCell->GetCellId() == 1)
             	 {
             		 //std::cout << "area:" << cell_target_area << "transition area " << 0.5+pCell->GetCellCycleModel()->GetG2Duration()/(cell_slope*4.0) << '\n' ;
             		 std::cout << "volume : " << cell_target_volume << "  G1 :  " << g1_duration << '\n' ;
             	 }
     } 
  
     // set cell data here 
     pCell->GetCellData()->SetItem("target volume", cell_target_volume); 
 	 
 } 
 	 
 template<unsigned DIM> 
 double TargetVolumeGrowthModifier<DIM>::GetMatureCellTargetVolume() 
 { 
     return mMatureCellTargetVolume; 
 } 
 	 
 template<unsigned DIM> 
 void TargetVolumeGrowthModifier<DIM>::SetMatureCellTargetVolume(double matureCellTargetVolume) 
 { 
     assert(matureCellTargetVolume >= 0.0); 
     mMatureCellTargetVolume = matureCellTargetVolume; 
 } 

 template<unsigned DIM>
 double TargetVolumeGrowthModifier<DIM>::GetMeanSlope()
 {
     return mMeanSlope;
 }

 template<unsigned DIM>
 void TargetVolumeGrowthModifier<DIM>::SetMeanSlope(double meanSlope)
 {
     assert(meanSlope >= 0.0);
     mMeanSlope = meanSlope;
 }
 ///////////////////////////////////////////////////////////////////////////// 
 // Explicit instantiation 
 ///////////////////////////////////////////////////////////////////////////// 
 	 
 template class TargetVolumeGrowthModifier<1>; 
 template class TargetVolumeGrowthModifier<2>; 
 template class TargetVolumeGrowthModifier<3>; 
 	 
 // Serialization for Boost >= 1.36 
 #include "SerializationExportWrapperForCpp.hpp" 
 EXPORT_TEMPLATE_CLASS_SAME_DIMS(TargetVolumeGrowthModifier)
