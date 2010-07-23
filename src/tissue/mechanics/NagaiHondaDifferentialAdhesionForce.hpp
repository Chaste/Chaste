/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#ifndef NAGAIHONDADIFFERENTIALADHESIONFORCE_HPP_
#define NAGAIHONDADIFFERENTIALADHESIONFORCE_HPP_


#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "VertexBasedTissue.hpp"

/**
 * The possible paired contact types types of TissueCells.
 */
typedef enum CellContactsType
{
    WILD_WILD,
    LABELLED_LABELLED,
    WILD_LABELLED,
    OTHER
} CellContactsType;

static const unsigned NUM_CELL_CONTACTS_TYPES=4;

/**
 * A force class for use in vertex-based tissue simulations, based on a mechanical
 * model proposed by T. Nagai and H. Honda ("A dynamic cell model for the formation
 * of epithelial tissues", Philosophical Magazine Part B 81:699-719).
 * Assume differential adhesion.
 */
template<unsigned DIM>
class NagaiHondaDifferentialAdhesionForce  : public AbstractForce<DIM>
{
friend class TestDifferentialAdhesionForcesNotForRelease;

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
    }

public:

    /**
     * Constructor.
     */
    NagaiHondaDifferentialAdhesionForce();

    /**
     * Destructor.
     */
    ~NagaiHondaDifferentialAdhesionForce();

    /**
     * Overridden AddForceContribution() method.
     *
     * Calculates the force on each node in the vertex-based tissue based on the
     * Nagai Honda model.
     *
     * @param rForces reference to vector of forces on nodes
     * @param rTissue reference to the tissue
     */
    void AddForceContribution(std::vector<c_vector<double, DIM> >& rForces, AbstractTissue<DIM>& rTissue);

    /**
     * Get the adhesion parameter for the edge between two given nodes.
     *
     * \todo We now are only testing 2 mutation states. This method could be extended
     * to handle any number of mutation states. One possibility is to have a method to set
     * a dictionary, where given a pair of cell types, the dictionary returns a corresponding
     * adhesion parameter value.
     *
     * \todo omit combinationCellType parameter and have method GetAdhesionParameterDifferentialAddition
     * call combinationCellType to get combinationCellType. Thus, make GetCombinationCellTypes private.
     *
     * @param pNodeA one node
     * @param pNodeB the other node
     * @param combinationCellType
     *
     * @return the adhesion parameter for this edge.
     */
    double GetAdhesionParameterDifferentialAddition(Node<DIM>* pNodeA, Node<DIM>* pNodeB, CellContactsType combinationCellType);

    /**
     * Get the combinationCellType
     *
     * \todo Make GetCombinationCellTypes private. See todo comment with GetAdhesionParameterDifferentialAddition.
     *
     * @param pNodeA one node
     * @param pNodeB the other node
     * @param rTissue reference to the tissue
     *
     * @return the combinationCellType for this edge.
     */
    CellContactsType GetCombinationCellTypes(Node<DIM>* pNodeA, Node<DIM>* pNodeB, AbstractTissue<DIM>& rTissue);
};


#include "SerializationExportWrapper.hpp"

EXPORT_TEMPLATE_CLASS_SAME_DIMS(NagaiHondaDifferentialAdhesionForce)

#endif /*NAGAIHONDADIFFERENTIALADHESIONFORCE_HPP_*/
