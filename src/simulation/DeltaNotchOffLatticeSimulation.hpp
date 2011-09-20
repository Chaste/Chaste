/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef DELTANOTCHOFFLATTICESIMULATION_HPP_
#define DELTANOTCHOFFLATTICESIMULATION_HPP_

#include <map>
#include "ChasteSerialization.hpp"
#include "OffLatticeSimulation.hpp"
#include "PetscTools.hpp"
#include "OutputFileHandler.hpp"

/**
 * Subclass of OffLatticeSimulation in which the mean levels Delta in neighbouring cells
 * are computed and stored in CellwiseData for use in DeltaNotchOdeSystem in a centre based cell population
 */
template<unsigned DIM>
class DeltaNotchOffLatticeSimulation : public OffLatticeSimulation<DIM>
{
private :

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * Serialization of singleton objects must be done with care.
     * Before the object is serialized via a pointer, it *MUST* be
     * serialized directly, or an assertion will trip when a second
     * instance of the class is created on de-serialization.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<OffLatticeSimulation<DIM> >(*this);
    }


    /**
     *  Overridden PostSolve() method.
     */
    void PostSolve();


public:

    /** The file that the values of beta catenin is written out to. */
	out_stream mVizDeltaFile;
    /**
     * Default constructor.
     *
     * @param rCellPopulation A cell population facade class (contains a mesh and cells)
     * @param deleteCellPopulationInDestructor whether to delete cell population and force collection.
     * @param initialiseCells whether to initialise cells (set to false when loading from an archive)
     */
     DeltaNotchOffLatticeSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
                                   bool deleteCellPopulationInDestructor=false,
                                   bool initialiseCells=true);

     /**
      * Destructor.
      */
    ~DeltaNotchOffLatticeSimulation();
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DeltaNotchOffLatticeSimulation)

namespace boost
{
namespace serialization
{
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const DeltaNotchOffLatticeSimulation<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    const AbstractCellPopulation<DIM> * p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, DeltaNotchOffLatticeSimulation<DIM> * t, const unsigned int file_version)
{
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    ::new(t)DeltaNotchOffLatticeSimulation<DIM>(*p_cell_population, true, false);
}
}
} // namespace ...

#endif /*DELTANOTCHOFFLATTICESIMULATION_HPP_*/
