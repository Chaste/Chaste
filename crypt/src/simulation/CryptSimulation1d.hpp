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

#ifndef CRYPTSIMULATION1D_HPP_
#define CRYPTSIMULATION1D_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include <cmath>
#include <iostream>

#include "WntConcentration.hpp"
#include "CryptSimulationBoundaryCondition.hpp"
#include "OffLatticeSimulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CryptCentreBasedDivisionRule.hpp"

/**
 * A 1D crypt simulation object. The model is a simplified version of a 2D crypt model
 * developed by Meineke et al (doi:10.1046/j.0960-7722.2001.00216.x).
 */
class CryptSimulation1d : public OffLatticeSimulation<1>
{
    // Allow tests to access private members to test private functions
    friend class TestCryptSimulation1d;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the simulation and member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<OffLatticeSimulation<1> >(*this);

        SerializableSingleton<WntConcentration<1> >* p_wnt_wrapper = WntConcentration<1>::Instance()->GetSerializationWrapper();
        archive & p_wnt_wrapper;
    }

    /** Helper member that is a static cast of the cell population. */
    MeshBasedCellPopulation<1>* mpStaticCastCellPopulation;

public:

    /**
     * Constructor.
     *
     * @param rCellPopulation A cell population object
     * @param deleteCellPopulationInDestructor Whether to delete the cell population on destruction to
     *     free up memory (defaults to false)
     * @param initialiseCells whether to initialise cells (defaults to true, set to false when loading from an archive)
     */
    CryptSimulation1d(AbstractCellPopulation<1>& rCellPopulation,
                      bool deleteCellPopulationInDestructor=false,
                      bool initialiseCells=true);

    /**
     * Destructor.
     *
     * This frees the CryptSimulationBoundaryCondition.
     */
    virtual ~CryptSimulation1d();

    /**
     * Outputs simulation parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CryptSimulation1d)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CryptSimulation1d.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const CryptSimulation1d * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<1>* p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise a CryptSimulation1d.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, CryptSimulation1d * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<1>* p_cell_population;
    ar & p_cell_population;

    // Invoke inplace constructor to initialise instance, last two variables set extra
    // member variables to be deleted as they are loaded from archive and to not initialise sells.
    ::new(t)CryptSimulation1d(*p_cell_population, true, false);
}
}
} // namespace

#endif /*CRYPTSIMULATION1D_HPP_*/
