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

#ifndef PARALLELCELLSGENERATOR_HPP_
#define PARALLELCELLSGENERATOR_HPP_

#include <vector>
#include "Cell.hpp"
#include "WildTypeCellMutationState.hpp"
#include "AbstractCellProperty.hpp"
#include "NodesOnlyMesh.hpp"

/**
 * A helper class for generating a vector of cells and a NodesOnlyMesh from an
 * archive file without loading all cells to memory on each processor.
 */
template<class CELL_CYCLE_MODEL, unsigned DIM>
class ParallelCellsGenerator
{
    friend class TestParallelCellsGenerator;    // For testing.

private:

    /**
     * Calculate the bounding box of all positions in the archive path.
     *
     * @param archivePath the path to the cell location archive.
     * @return the bounding box in format [xmin, xmax...].
     */
    c_vector<double, 2*DIM> GetArchiveBoundingBox(std::string archivePath);

public:

    /**
     * Generate nodes and cells for a given decomposition of space, dictated by mesh.
     * @param archivePath the name of the file containing the cell locations
     * @param cells a reference to an empty vector into which cells can be inserted.
     * @param mesh a reference to the mesh that will store the nodes.
     * @param pCellProliferativeType the proliferative type of the cells.
     */
    void GenerateParallelCells( std::string archivePath,
                                std::vector<CellPtr>& cells,
                                NodesOnlyMesh<DIM>& mesh,
                                boost::shared_ptr<AbstractCellProperty> pCellProliferativeType=boost::shared_ptr<AbstractCellProperty>());
};

// Implementation
template<class CELL_CYCLE_MODEL, unsigned DIM>
void ParallelCellsGenerator<CELL_CYCLE_MODEL, DIM>::GenerateParallelCells(std::string archivePath, std::vector<CellPtr>& cells, NodesOnlyMesh<DIM>& mesh, boost::shared_ptr<AbstractCellProperty> pCellProliferativeType)
{
    // Get a bounding box for the archived nodes
    c_vector<double, 2*DIM> base_bounding_box = GetArchiveBoundingBox(archivePath);

    // Add small factor to make sure the box is big enough to contain the initial population.
    for (unsigned i=0; i<DIM; i++)
    {
        base_bounding_box[2*i] -= 1e-14;
        base_bounding_box[2*i + 1] += 1e-14;
    }

    // Construct box collection of NodesOnlyMesh
    mesh.SetInitialBoxCollection(base_bounding_box, mesh.GetMaximumInteractionDistance());

    unsigned node_index = 0;

    std::ifstream infile(archivePath.c_str(), std::ios::in);

    std::string line;
    std::getline(infile, line);    // Get first line which is ignored as it is a header.

    while (std::getline(infile, line))
    {
        // Get line in archive file.
        std::istringstream iss(line);

        // Extract a location.
        c_vector<double, DIM> location;
        for (unsigned k=0; k<DIM; k++)
        {
            iss >> location[k];
        }

        // Add it to the mesh and make an appropriate cell if it is owned.
        if (mesh.IsOwned(location))
        {
            // Make a node at this point.
            Node<DIM>* p_node = new Node<DIM>(node_index, location);    // Memory is managed by NodesOnlyMesh.
            mesh.AddNode(p_node);

            // Make cell
            CELL_CYCLE_MODEL* p_cell_cycle_model = new CELL_CYCLE_MODEL;
            p_cell_cycle_model->SetDimension(DIM);

            // Wild type hard-coded for now.
            boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
            p_cell->SetCellProliferativeType(pCellProliferativeType);
            cells.push_back(p_cell);

            node_index++;
        }
    }

    infile.close();

    // Tidy up the mesh.
    NodeMap map(mesh.GetMaximumNodeIndex()+1);
    mesh.ReMesh(map);
}

template<class CELL_CYCLE_MODEL, unsigned DIM>
c_vector<double, 2*DIM> ParallelCellsGenerator<CELL_CYCLE_MODEL, DIM>::GetArchiveBoundingBox(std::string archivePath)
{
    c_vector<double, 2*DIM> bounding_box;
    for (unsigned i=0; i<DIM; i++)
    {
        bounding_box[2*i] = DBL_MAX;
        bounding_box[2*i + 1] = -DBL_MAX;
    }

    std::ifstream infile(archivePath.c_str());

    std::string line;
    std::getline(infile, line);

    std::istringstream iss(line);

    // Check file dimensions and template dimensions match.
    unsigned num_lines;
    unsigned file_dimension;

    iss >> num_lines >> file_dimension;

    if (file_dimension != DIM)
    {
        EXCEPTION("Space dimension of ParallelCellsGenerator and archive file do not match");
    }

    while (std::getline(infile, line))
    {
        std::stringstream new_iss(line);
        for (unsigned i=0; i<DIM; i++)
        {
            double point;
            new_iss >> point;
            bounding_box[2*i] = (point < bounding_box[2*i]) ? point : bounding_box[2*i];
            bounding_box[2*i+1] = (point > bounding_box[2*i+1]) ? point : bounding_box[2*i+1];
        }
    }

    infile.close();

    return bounding_box;
}
#endif /* PARALLELCELLSGENERATOR_HPP_ */
