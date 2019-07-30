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

#ifndef ABSTRACTPURKINJECELLFACTORY_HPP_
#define ABSTRACTPURKINJECELLFACTORY_HPP_

#include "AbstractCardiacCellFactory.hpp"
#include "MixedDimensionMesh.hpp"
#include "FakeBathCell.hpp"
/**
 *  Subclass for also creating cell models for Purkinje cells.
 *  The user has to implement the CreatePurkinjeCellForTissueNode() method.
 *
 *  The dimensions should be 2 or 3.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class AbstractPurkinjeCellFactory : public AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>
{
private:
    /**
     * Reads in node id and resistance values of junction nodes
     *
     * The .pvj file path is specified by HeartConfig::Instance()->GetMeshFileName() + ".pvj"
     *
     * Note, this is called by SetMesh
     */
    void ReadJunctionsFile();

protected:
    /** Saved pointer to the mixed dimension mesh */
    MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>* mpMixedDimensionMesh;

    /** A set of local purkinje node indices */
    std::set<unsigned> mLocalPurkinjeNodes;

    /** A map between junction node ids and resistances */
    std::map<unsigned, double> mJunctionMap;

    /**
     * @return a newly created purkinje cells for the given node.
     * Must be overridden by subclasses to return a Purkinje cell object for the given node.
     * @param pNode  pointer to node object
     * @param pCardiacCell  the cardiac cell that has already been created at this node
     */
    virtual AbstractCardiacCellInterface* CreatePurkinjeCellForTissueNode(Node<SPACE_DIM>* pNode,
                                                                          AbstractCardiacCellInterface* pCardiacCell)=0;

    /**
     * @return a newly created purkinje-ventricular junction between the two cells provided.
     *
     * @param pNode  the node in the mesh at which this junction is located
     * @param pPurkinjeCell  the Purkinje cell
     * @param pCardiacCell  the ventricular cell
     * @param resistance  the junction resistance, in kilo-Ohms
     */
    void CreateJunction(const Node<SPACE_DIM>* pNode,
                        AbstractCardiacCellInterface* pPurkinjeCell,
                        AbstractCardiacCellInterface* pCardiacCell,
                        double resistance);


    /**
     * @return a newly created purkinje-ventricular junction between the two cells provided if the junction is defined in the corresponding .pvj file.
     *
     * ReadJunctionsFile() must be called before calling this method
     *
     * @param pNode  the node in the mesh at which this junction is located
     * @param pPurkinjeCell  the Purkinje cell
     * @param pCardiacCell  the ventricular cell
     */
    void CreateJunctionFromFile(const Node<SPACE_DIM>* pNode,
                                AbstractCardiacCellInterface* pPurkinjeCell,
                                AbstractCardiacCellInterface* pCardiacCell);


public:

    /* Constructor does nothing */
    AbstractPurkinjeCellFactory();


    /** Overridden set mesh which must take a MixedDimensionMesh
     *  @param pMesh Pointer to the mesh. */
    void SetMesh(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh);

    /**
     * @return a newly created cell object for the given node.
     *
     * The default implementation checks whether the node is in a Purkinje node, in which
     * case it calls CreatePurkinjeCellForTissueNode (which must be defined by subclasses),
     * otherwise it returns a pointer to a (unique) fake cell
     *
     * @param pNode  pointer to node
     * @param pCardiacCell  the cardiac cell that has already been created at this node
     */
    AbstractCardiacCellInterface* CreatePurkinjeCellForNode(Node<SPACE_DIM>* pNode,
                                                            AbstractCardiacCellInterface* pCardiacCell);

    /**
     * May be overridden by subclasses to perform any necessary work after all Purkinje cells
     * have been created.
     *
     * @param pPurkinjeCellsDistributed  Pointer to a vector of Purkinje cell pointers.
     * @param lo  Lowest index owned by this process.
     * @param hi  Highest index owned by this process.
     */
    virtual void FinalisePurkinjeCellCreation(std::vector< AbstractCardiacCellInterface* >* pPurkinjeCellsDistributed,
                                              unsigned lo, unsigned hi)
    {
    }

    /**
     *  @return the mixed dimension mesh (for possible use in CreatePurkinjeCellForTissueNode()).
     *  Note: GetMesh() just returns a pointer to an AbstractTetrahedralMesh.
     */
    MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>* GetMixedDimensionMesh();
};

#endif // ABSTRACTPURKINJECELLFACTORY_HPP_
