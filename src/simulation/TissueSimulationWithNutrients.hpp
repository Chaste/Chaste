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
#ifndef TISSUESIMULATIONWITHNUTRIENTS_HPP_
#define TISSUESIMULATIONWITHNUTRIENTS_HPP_

#include <map>
#include "ChasteSerialization.hpp"

#include "TissueSimulation.hpp"

#include "AbstractLinearEllipticPde.hpp"
#include "AveragedSinksPde.hpp"
#include "TetrahedralMesh.hpp"
#include "PetscTools.hpp"

#include "OutputFileHandler.hpp"
#include "TissueCell.hpp"

/**
 * A nutrient-dependent tissue simulation class.
 */
template<unsigned DIM>
class TissueSimulationWithNutrients : public TissueSimulation<DIM>
{
    // Allow tests to access private members, in order to test computation of private functions
    friend class TestTissueSimulationWithNutrients;

private :

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & boost::serialization::base_object<TissueSimulation<DIM> >(*this);
        archive & mWriteAverageRadialNutrientResults;
        archive & mWriteDailyAverageRadialNutrientResults;
        archive & mNumRadialIntervals;
        archive & mCellNutrientElementMap;
    }

    /**
     *  Current nutrient concentration, for use as an initial guess
     *  when solving the nutrient PDE.
     */
    Vec mNutrientSolution;

    /**
     *  Pointer to the PDE satisfied by the nutrient.
     */
    AbstractLinearEllipticPde<DIM,DIM>* mpPde;

    /**
     *  Pointer to the averaged sink PDE satisfied by the nutrient.
     */
    AveragedSinksPde<DIM>* mpAveragedSinksPde;

    /**
     *  File that the nutrient values are written out to.
     */
    out_stream mpNutrientResultsFile;

    /**
     *  File that the average radial nutrient distribution is written out to.
     */
    out_stream mpAverageRadialNutrientResultsFile;

    /**
     *  Whether to write to file the average radial nutrient distribution.
     */
    bool mWriteAverageRadialNutrientResults;

    /**
     *  Whether to write the average radial nutrient distribution DAILY.
     */
    bool mWriteDailyAverageRadialNutrientResults;

    /**
     *  Number of radial 'bins' used to calculate the average
     *  radial nutrient distribution.
     */
    unsigned mNumRadialIntervals;

    /**
     *  Coarse nutrient mesh on which to solve the nutrient PDE.
     */
    TetrahedralMesh<DIM,DIM>* mpCoarseNutrientMesh;

    /**
     * Map between cells and the elements of the coarse nutrient mesh containing them.
     */
    std::map<TissueCell*, unsigned> mCellNutrientElementMap;

    /**
     *  Overridden SetupSolve() method.
     */
    void SetupSolve();

    /**
     *  Set up the nutrient writer.
     */
    void SetupWriteNutrient();

    /**
     * Write the nutrient distribution to file at a specified time.
     *
     * @param time The time at which to record the nutrient distribution
     */
    void WriteNutrient(double time);

    /**
     * Write the average radial nutrient distribution to file at a specified time.
     *
     * @param time The time at which to record the average radial nutrient distribution
     * @param numIntervals  The number of radial intervals in which the average nutrient concentration is calculated
     */
    void WriteAverageRadialNutrientDistribution(double time, unsigned numIntervals);

    /**
     *  Solve the nutrient PDE.
     */
    void SolveNutrientPde();

    /**
     *  Solve the nutrient PDE on a coarse mesh.
     */
    void SolveNutrientPdeUsingCoarseMesh();

    /**
     * Find the index of the coarse mesh element containing a given cell.
     *
     * @param rCell the cell
     *
     * @return the element index.
     */
    unsigned FindElementContainingCell(TissueCell& rCell);

    /**
     *  Overridden PostSolve() method.
     */
    void PostSolve();

    /**
     *  Overridden AfterSolve() method.
     */
    void AfterSolve();

    /**
     *  Create a coarse mesh on which to solve the nutrient PDE.
     *
     * \todo currently only works in 2D (see #737)
     *
     * @param coarseGrainScaleFactor the ratio of the width of the coarse nutrient mesh to the initial width of the tissue
     */
    void CreateCoarseNutrientMesh(double coarseGrainScaleFactor);

    /**
     *  Initialise the std::map mCellNutrientElementMap.
     */
    void InitialiseCoarseNutrientMesh();

    /**
     * Overridden WriteVisualizerSetupFile() method.
     *
     * Writes out special information about the mesh to the visualizer.
     */
    void WriteVisualizerSetupFile();

public:

    /**
     * Constructor
     *
     * @param rTissue A tissue facade class (contains a mesh and cells)
     * @param forceCollection The mechanics to use in the simulation
     * @param pPde The PDE for the nutrient concentration(s)
     * @param pAveragedSinksPde The PDE for the nutrient concentration(s)
     * @param deleteTissueAndForceCollection whether to delete the tissue on destruction to free up memory
     * @param initialiseCells whether to initialise cells (set to false when loading from an archive)
     *
     */
     TissueSimulationWithNutrients(AbstractTissue<DIM>& rTissue,
                                   std::vector<AbstractForce<DIM>*> forceCollection,
                                   AbstractLinearEllipticPde<DIM,DIM>* pPde=NULL,
                                   AveragedSinksPde<DIM>* pAveragedSinksPde=NULL,
                                   bool deleteTissueAndForceCollection=false,
                                   bool initialiseCells=true);

    /**
     * Destructor
     *
     * Free any memory allocated by the constructor.
     * This frees the current nutrient distribution, if it exists.
     */
    ~TissueSimulationWithNutrients();

    /**
     * A small hack until we fully archive this class -
     * needed to set the PDE after loading a simulation
     * from an archive.
     *
     * @param pPde pointer to the PDE object
     */
    void SetPde(AbstractLinearEllipticPde<DIM,DIM>* pPde);

    /**
     * A small hack until we fully archive this class -
     * needed to set the PDE after loading a simulation
     * from an archive.
     *
     * @param pAveragedSinksPde pointer to the PDE object
     */
    void SetAveragedSinksPde(AveragedSinksPde<DIM>* pAveragedSinksPde);

    /**
     *  Get the current nutrient solution
     */
    Vec GetNutrientSolution();

    /**
     * Write the final (and optionally also the daily) average
     * radial nutrient distribution to file.
     *
     * @param numRadialIntervals The number of radial intervals in which the average nutrient concentration is calculated
     * @param writeDailyResults Whether to record the average radial nutrient distribution at the end of each day of the simulation
     */

    void SetWriteAverageRadialNutrientResults(unsigned numRadialIntervals=10,
                                              bool writeDailyResults=false);

    /**
     * Solve the nutrient PDE on a coarse mesh.
     *
     * @param coarseGrainScaleFactor the ratio of the width of the coarse nutrient mesh to the initial width of the tissue
     */
    void UseCoarseNutrientMesh(double coarseGrainScaleFactor=10.0);

};


#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(TissueSimulationWithNutrients)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a TissueSimulationWithNutrients.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const TissueSimulationWithNutrients<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractTissue<DIM> * p_tissue = &(t->rGetTissue());
    ar & p_tissue;
    const std::vector<AbstractForce<DIM>*> force_collection = t->rGetForceCollection();
    ar & force_collection;
}

/**
 * De-serialize constructor parameters and initialise a TissueSimulationWithNutrients.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, TissueSimulationWithNutrients<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractTissue<DIM>* p_tissue;
    ar >> p_tissue;
    std::vector<AbstractForce<DIM>*> force_collection;
    ar >> force_collection;

    // Invoke inplace constructor to initialise instance
    ::new(t)TissueSimulationWithNutrients<DIM>(*p_tissue, force_collection, NULL, NULL, true, false);
}
}
} // namespace ...


#endif /*TISSUESIMULATIONWITHNUTRIENTS_HPP_*/
