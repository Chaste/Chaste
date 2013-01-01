/*

Copyright (c) 2005-2013, University of Oxford.
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

#ifndef VOLUMEDEPENDENTAVERAGEDSOURCEPDE_HPP_
#define VOLUMEDEPENDENTAVERAGEDSOURCEPDE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "NodeBasedCellPopulation.hpp"
#include "AveragedSourcePde.hpp"
#include "TetrahedralMesh.hpp"
#include "AbstractLinearEllipticPde.hpp"

/**
 *  A PDE which calculates the source term by adding the number of cells
 *  in the element containing that point and scaling by the element area.
 */
template<unsigned DIM>
class VolumeDependentAveragedSourcePde : public AveragedSourcePde<DIM>
{
    friend class TestCellBasedPdes;

private:

    /** Needed for serialization.*/
    friend class boost::serialization::access;
    /**
     * Serialize the PDE and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
       archive & boost::serialization::base_object<AveragedSourcePde<DIM> >(*this);
    }

    /** Static cast of the NodeBasedCellPopulation. **/
    NodeBasedCellPopulation<DIM>* mpStaticCastCellPopulation;

public:

    /**
     * Constructor.
     *
     * @param rCellPopulation reference to the cell population
     * @param coefficient the coefficient of consumption of nutrient by cells (defaults to 0.0)
     */
    VolumeDependentAveragedSourcePde(AbstractCellPopulation<DIM>& rCellPopulation, double coefficient=0.0);

    /**
     * Set up the source terms.
     *
     * @param rCoarseMesh reference to the coarse mesh
     * @param pCellPdeElementMap optional pointer to the map from cells to coarse elements
     */
    void SetupSourceTerms(TetrahedralMesh<DIM,DIM>& rCoarseMesh, std::map<CellPtr, unsigned>* pCellPdeElementMap=NULL);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VolumeDependentAveragedSourcePde)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a VolumeDependentAveragedSourcePde.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const VolumeDependentAveragedSourcePde<DIM>* t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise a VolumeDependentAveragedSourcePde.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, VolumeDependentAveragedSourcePde<DIM>* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)VolumeDependentAveragedSourcePde<DIM>(*p_cell_population);
}
}
} // namespace ...

#endif /*VOLUMEDEPENDENTAVERAGEDSOURCEPDE_HPP_*/
