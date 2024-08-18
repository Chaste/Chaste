/*

Copyright (c) 2005-2024, University of Oxford.
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

/*

Copyright (c) 2005-2024, University of Oxford.
 All rights reserved.

 University of Oxford means the Chancellor, Masters and Scholars of the
 University of Oxford, having an administrative office at Wellington
 Square, Oxford OX1 2JD, UK.

 This file is CellPopulation of Chaste.

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
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A CellPopulationICULAR PURPOSE
 ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#ifndef CELLPOPULATIONPYCHASTEACTORGENERATOR_HPP_
#define CELLPOPULATIONPYCHASTEACTORGENERATOR_HPP_

#include <vector>
#include "SmartPointers.hpp"
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkLookupTable.h>
#include "AbstractCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "AbstractPyChasteActorGenerator.hpp"

/**
 * This class generates VTK actors for CellPopulations
 */
template<unsigned DIM>
class CellPopulationPyChasteActorGenerator : public AbstractPyChasteActorGenerator<DIM>
{
    /**
     * The CellPopulation
     */
    boost::shared_ptr<AbstractCellPopulation<DIM> > mpCellPopulation;

    /**
     * Show mutable mesh edges for Mesh Based populations
     */
    bool mShowMutableMeshEdges;

    /**
     * Show voronoi mesh edges for Mesh and Vertex Based populations
     */
    bool mShowVoronoiMeshEdges;

    /**
     * Show Potts mesh edges for Ca and Potts Based populations
     */
    bool mShowPottsMeshEdges;

    /**
     * Show Potts mesh outlines
     */
    bool mShowPottsMeshOutlines;

    /**
     * Color the cells by type
     */
    bool mColorByCellType;

    /**
     * Color the cells by data
     */
    bool mColorByCellData;

    /**
     * Color the cells by mutation state
     */
    bool mColorByCellMutationState;

    /**
     * Color the cells by label
     */
    bool mColorByCellLabel;

    /**
     * Whether to show the cell centres
     */
    bool mShowCellCentres;

    /**
     * Color cells using a user defined color
     */
    bool mColorCellByUserDefined;

public:

    /**
     * Constructor
     */
    CellPopulationPyChasteActorGenerator();

    /**
     * Destructor
     */
    ~CellPopulationPyChasteActorGenerator();

    /**
     * Add the CellPopulation actor to the renderer
     * @param pRenderer the current renderer
     */
    void AddActor(vtkSmartPointer<vtkRenderer> pRenderer);

    /**
     * Specialized class for adding Mesh based population
     * @param pRenderer the current renderer
     */
    void AddMeshBasedCellPopulationActor(vtkSmartPointer<vtkRenderer> pRenderer);

    /**
     * Specialized class for adding Mesh based population
     * @param pRenderer the current renderer
     */
    void AddVertexBasedCellPopulationActor(vtkSmartPointer<vtkRenderer> pRenderer);

    /**
     * Specialized class for adding Immersed boundary population
     * @param pRenderer the current renderer
     */
    void AddImmersedBoundaryCellPopulationActor(vtkSmartPointer<vtkRenderer> pRenderer);

    /**
     * Specialized class for adding Ca based population
     * @param pRenderer the current renderer
     */
    void AddCaBasedCellPopulationActor(vtkSmartPointer<vtkRenderer> pRenderer);

    /**
     * Specialized class for adding Potts based population
     * @param pRenderer the current renderer
     */
    void AddPottsBasedCellPopulationActor(vtkSmartPointer<vtkRenderer> pRenderer);

    /**
     * Set the CellPopulation
     * @param pCellPopulation the CellPopulation to render
     */
    void SetCellPopulation(boost::shared_ptr<AbstractCellPopulation<DIM> > pCellPopulation);

    /**
     * @param showEdges show the voronoi mesh
     */
    void SetShowVoronoiMeshEdges(bool showEdges);

    /**
     * @param showEdges show the mutable mesh
     */
    void SetShowMutableMeshEdges(bool showEdges);

    /**
     * @param showEdges show the potts mesh
     */
    void SetShowPottsMeshEdges(bool showEdges);

    /**
     * @param showOutlines show the outlines of Potts cells
     */
    void SetShowPottsMeshOutlines(bool showOutlines);

    /**
     * @param colorByCellType color cells by type
     */
    void SetColorByCellType(bool colorByCellType);

    /**
     * @param colorByCellType color cells by mutation state
     */
    void SetColorByCellMutationState(bool colorByCellMutationState);

    /**
     * @param colorByCellType color cells by label
     */
    void SetColorByCellLabel(bool colorByCellLabel);

    /**
     * @param colorByCellType color cells by the user defined point color
     */
    void SetColorByUserDefined(bool colorByCellUserDefined);

    /**
     * @param colorByCellData color cells by data
     */
    void SetColorByCellData(bool colorByCellData);

    /**
     * @param showCentres show cell centres
     */
    void SetShowCellCentres(bool showCentres);
};

#endif /* CELLPOPULATIONPYCHASTEACTORGENERATOR_HPP_*/
