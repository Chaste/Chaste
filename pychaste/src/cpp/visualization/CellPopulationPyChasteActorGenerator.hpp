/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef CELLPOPULATIONPYCHASTEACTORGENERATOR_HPP_
#define CELLPOPULATIONPYCHASTEACTORGENERATOR_HPP_

#include <vector>

#include <vtkLookupTable.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>

#include "AbstractCellPopulation.hpp"
#include "AbstractPyChasteActorGenerator.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "SmartPointers.hpp"
#include "VertexBasedCellPopulation.hpp"

class vtkDataArray;
class vtkPoints;
class vtkUnstructuredGrid;

/**
 * Generates VTK actors that render a Chaste cell population for visualisation. 
 *
 * Given a population set via SetCellPopulation(), AddActor() builds the actors
 * appropriate to the population type (mesh-, vertex-, Potts-, CA-, node- or
 * immersed-boundary-based) and adds them to a renderer.
 * What is drawn, how cells are coloured, and which mesh edges or cell centres
 * to show are configured through the Set... methods and the options inherited
 * from AbstractPyChasteActorGenerator.
 */
template <unsigned DIM>
class CellPopulationPyChasteActorGenerator : public AbstractPyChasteActorGenerator<DIM>
{
    /**
     * The CellPopulation
     */
    boost::shared_ptr<AbstractCellPopulation<DIM> > mpCellPopulation;

    /**
     * Show mutable mesh edges for Mesh Based populations
     */
    bool mShowMutableMeshEdges = false;

    /**
     * Show voronoi mesh edges for Mesh and Vertex Based populations
     */
    bool mShowVoronoiMeshEdges = true;

    /**
     * Show Potts mesh edges for Ca and Potts Based populations
     */
    bool mShowPottsMeshEdges = false;

    /**
     * Show Potts mesh outlines
     */
    bool mShowPottsMeshOutlines = false;

    /**
     * Color the cells by type
     */
    bool mColorByCellType = false;

    /**
     * Color the cells by data
     */
    bool mColorByCellData = false;

    /**
     * Color the cells by mutation state
     */
    bool mColorByCellMutationState = false;

    /**
     * Color the cells by label
     */
    bool mColorByCellLabel = false;

    /**
     * Whether to show the cell centres
     */
    bool mShowCellCentres = false;

    /**
     * Color cells using a user defined color
     */
    bool mColorCellByUserDefined = false;

public:

    /**
     * Constructor
     */
    CellPopulationPyChasteActorGenerator() = default;

    /**
     * Destructor
     */
    ~CellPopulationPyChasteActorGenerator() = default;

    /**
     * Build the actors for the current cell population and add them to the
     * renderer. Optionally draws cell centres, then dispatches on the population
     * type to the matching helper below.
     *
     * @param pRenderer the renderer to add the actors to
     */
    void AddActor(vtkSmartPointer<vtkRenderer> pRenderer);

    /**
     * Add the actors for a CA-based cell population to the renderer.
     * Normally called via AddActor() once the population type has been determined.
     *
     * @param pRenderer the renderer to add the actors to
     */
    void AddCaBasedCellPopulationActor(vtkSmartPointer<vtkRenderer> pRenderer);

    /**
     * Add the actors for an immersed-boundary cell population to the renderer.
     * Normally called via AddActor() once the population type has been determined.
     *
     * @param pRenderer the renderer to add the actors to
     */
    void AddImmersedBoundaryCellPopulationActor(vtkSmartPointer<vtkRenderer> pRenderer);

    /**
     * Add the actors for a mesh-based cell population to the renderer: the
     * Voronoi tessellation and/or the underlying (Delaunay) mesh edges.
     * Normally called via AddActor() once the population type has been determined.
     *
     * @param pRenderer the renderer to add the actors to
     */
    void AddMeshBasedCellPopulationActor(vtkSmartPointer<vtkRenderer> pRenderer);

    /**
     * Add the actors for a Potts-based cell population to the renderer.
     * Normally called via AddActor() once the population type has been determined.
     *
     * @param pRenderer the renderer to add the actors to
     */
    void AddPottsBasedCellPopulationActor(vtkSmartPointer<vtkRenderer> pRenderer);

    /**
     * Add the actors for a vertex-based cell population to the renderer.
     * Normally called via AddActor() once the population type has been determined.
     *
     * @param pRenderer the renderer to add the actors to
     */
    void AddVertexBasedCellPopulationActor(vtkSmartPointer<vtkRenderer> pRenderer);

    /**
     * Set the cell population to render.
     * @param pCellPopulation the cell population to render
     */
    void SetCellPopulation(boost::shared_ptr<AbstractCellPopulation<DIM> > pCellPopulation);

    /**
     * Set whether to colour cells by a named CellData entry (see the data label
     * option inherited from AbstractPyChasteActorGenerator).
     * @param colorByCellData whether to colour cells by data
     */
    void SetColorByCellData(bool colorByCellData);

    /**
     * Set whether to colour cells by their CellLabel.
     * @param colorByCellLabel whether to colour cells by label
     */
    void SetColorByCellLabel(bool colorByCellLabel);

    /**
     * Set whether to colour cells by their mutation state.
     * @param colorByCellMutationState whether to colour cells by mutation state
     */
    void SetColorByCellMutationState(bool colorByCellMutationState);

    /**
     * Set whether to colour cells by their proliferative type.
     * @param colorByCellType whether to colour cells by type
     */
    void SetColorByCellType(bool colorByCellType);

    /**
     * Set whether to colour every cell with a single user-defined colour.
     * @param colorByCellUserDefined whether to use the user-defined colour
     */
    void SetColorByUserDefined(bool colorByCellUserDefined);

    /**
     * Set whether to draw cell centres as glyphs.
     * @param showCentres whether to show cell centres
     */
    void SetShowCellCentres(bool showCentres);

    /**
     * Set whether to draw the underlying (Delaunay) mesh edges, for mesh-based
     * populations.
     * @param showEdges whether to show the mutable mesh edges
     */
    void SetShowMutableMeshEdges(bool showEdges);

    /**
     * Set whether to draw the Potts/CA lattice edges, for Potts- and CA-based
     * populations.
     * @param showEdges whether to show the Potts mesh edges
     */
    void SetShowPottsMeshEdges(bool showEdges);

    /**
     * Set whether to draw outlines around Potts cells.
     * @param showOutlines whether to show the Potts cell outlines
     */
    void SetShowPottsMeshOutlines(bool showOutlines);

    /**
     * Set whether to draw the Voronoi tessellation edges, for mesh- and
     * vertex-based populations.
     * @param showEdges whether to show the Voronoi mesh edges
     */
    void SetShowVoronoiMeshEdges(bool showEdges);

private:

    /**
     * Render a coloured cell grid: adds a colour-mapped surface actor and
     * tube-filtered boundary edges to the renderer, plus a scale bar when
     * colouring by data. Shared by the vertex-, immersed-boundary- and
     * mesh-based populations.
     *
     * @param pRenderer the renderer to add actors to
     * @param pGrid the grid, carrying a "CellColors" cell-data array
     * @param thresholdGhosts whether to first threshold out negatively-tagged
     *        (ghost) cells, as used by mesh-based populations
     */
    void AddColouredVoronoiGrid(vtkSmartPointer<vtkRenderer> pRenderer,
                                vtkSmartPointer<vtkUnstructuredGrid> pGrid,
                                bool thresholdGhosts);

    /**
     * Append a location to a vtkPoints array, padding unused coordinates with
     * zero so that 1D and 2D locations are stored as 3D points.
     *
     * @param pPoints the points array to append to
     * @param rLocation the location to append
     */
    void AppendPoint(vtkPoints* pPoints, const c_vector<double, DIM>& rLocation) const;

    /**
     * Build a coloured VTK grid from the elements of a population's mesh, tagging
     * each cell with its colour value in a "CellColors" cell-data array. Shared by
     * the vertex-based and immersed-boundary populations, whose meshes expose the
     * same element-iteration interface.
     *
     * @param pPopulation the cell population
     * @return the coloured unstructured grid
     */
    template <typename PopulationType>
    vtkSmartPointer<vtkUnstructuredGrid> BuildMeshElementGrid(boost::shared_ptr<PopulationType> pPopulation);

    /**
     * Build a colour transfer function scaled to the range of a data array, using
     * the continuous (colour-by-data) or discrete colour transfer function.
     *
     * @param pArray the data array whose range the lookup is scaled to
     * @return the scaled colour transfer function
     */
    vtkSmartPointer<vtkColorTransferFunction> BuildScaledColourLookup(vtkDataArray* pArray) const;

    /**
     * Compute the scalar value used to colour a cell, according to the active
     * colour-by mode (type, data, mutation state, label, or cell id by default).
     *
     * @param pCell the cell to colour
     * @return the scalar colour value for the cell
     */
    double GetCellColour(const CellPtr& pCell) const;
};

#endif // CELLPOPULATIONPYCHASTEACTORGENERATOR_HPP_
