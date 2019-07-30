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

#ifndef VOLUMEDEPENDENTAVERAGEDSOURCEELLIPTICPDE_HPP_
#define VOLUMEDEPENDENTAVERAGEDSOURCEELLIPTICPDE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "NodeBasedCellPopulation.hpp"
#include "AveragedSourceEllipticPde.hpp"
#include "TetrahedralMesh.hpp"
#include "AbstractLinearEllipticPde.hpp"

/**
 * An elliptic PDE to be solved numerically using the finite element method, for
 * coupling to a cell-based simulation.
 *
 * This class inherits from AveragedSourceEllipticPde and may only be used with
 * a NodeBasedCellPopulation, since it assumes that each cell is associated with
 * a Node object.
 *
 * The PDE takes the form
 *
 * Grad.(D*Grad(u)) + k*u*rho(x) = 0,
 *
 * where the scalars D and k are specified by the members mDiffusionCoefficient and
 * mSourceCoefficient, respectively. Their values must be set in the constructor.
 *
 * The function rho(x) denotes the local density of non-apoptotic cells. This
 * quantity is computed for each element of a 'coarse' finite element mesh that is
 * passed to the method SetupSourceTerms() and stored in the member mCellDensityOnCoarseElements.
 *
 * For a point x, rho(x) is a weighted sum of the non-apoptotic cells whose centres
 * lie in each finite element containing that point, scaled by the area of that element.
 * The weighting assigned to each cell is given by the square of the radius of the
 * associated node, which is accessed using the GetRadius() method.
 *
 * \todo Consider creating a VolumeDependentAveragedSourceParabolicPde class
 */
template<unsigned DIM>
class VolumeDependentAveragedSourceEllipticPde : public AveragedSourceEllipticPde<DIM>
{
    friend class TestCellBasedEllipticPdes;

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
       archive & boost::serialization::base_object<AveragedSourceEllipticPde<DIM> >(*this);
    }

    /** Static cast of the NodeBasedCellPopulation. */
    NodeBasedCellPopulation<DIM>* mpStaticCastCellPopulation;

public:

    /**
     * Constructor.
     *
     * @param rCellPopulation reference to the cell population
     * @param coefficient the coefficient of consumption of nutrient by cells (defaults to 0.0)
     */
    VolumeDependentAveragedSourceEllipticPde(AbstractCellPopulation<DIM>& rCellPopulation, double coefficient=0.0);

    /**
     * Set up the source terms.
     *
     * @param rCoarseMesh reference to the coarse mesh
     * @param pCellPdeElementMap optional pointer to the map from cells to coarse elements
     */
    void SetupSourceTerms(TetrahedralMesh<DIM,DIM>& rCoarseMesh, std::map<CellPtr, unsigned>* pCellPdeElementMap=nullptr);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VolumeDependentAveragedSourceEllipticPde)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a VolumeDependentAveragedSourceEllipticPde.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const VolumeDependentAveragedSourceEllipticPde<DIM>* t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise a VolumeDependentAveragedSourceEllipticPde.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, VolumeDependentAveragedSourceEllipticPde<DIM>* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)VolumeDependentAveragedSourceEllipticPde<DIM>(*p_cell_population);
}
}
} // namespace ...

#endif /*VOLUMEDEPENDENTAVERAGEDSOURCEELLIPTICPDE_HPP_*/
