/*

Copyright (c) 2005-2017, University of Oxford.
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

#ifndef CELLWISESOURCEPARABOLICPDE_HPP_
#define CELLWISESOURCEPARABOLICPDE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellPopulation.hpp"
#include "AbstractLinearParabolicPdeSystem.hpp"

/**
 * A parabolic PDE to be solved numerically using the finite element method, for
 * coupling to a cell-based simulation.
 *
 * The PDE takes the form
 *
 * c*du/dt = Grad.(D*Grad(u)) + k*u*rho(x),
 *
 * where the scalars c, D and k are specified by the members mDuDtCoefficient,
 * mDiffusionCoefficient and mSourceCoefficient, respectively. Their values must
 * be set in the constructor.
 *
 * For a node of the finite element mesh with location x, the function rho(x)
 * equals one if there is a non-apoptotic cell associated with x, and
 * zero otherwise. Here, 'associated with' takes a different meaning for each
 * cell population class, and is encoded in the method IsPdeNodeAssociatedWithNonApoptoticCell().
 */
template<unsigned DIM>
class CellwiseSourceParabolicPde : public AbstractLinearParabolicPdeSystem<DIM, DIM>
{
    friend class TestCellBasedParabolicPdes;

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
       archive & boost::serialization::base_object<AbstractLinearParabolicPdeSystem<DIM, DIM> >(*this);
       archive & mDuDtCoefficient;
       archive & mDiffusionCoefficient;
       archive & mSourceCoefficient;
    }

protected:

    /** The cell population member. */
    AbstractCellPopulation<DIM, DIM>& mrCellPopulation;

    /** Coefficient of rate of change term.  */
    double mDuDtCoefficient;

    /** Diffusion coefficient. */
    double mDiffusionCoefficient;

    /** Coefficient of the rate of uptake of the dependent variable by non-apoptotic cells. */
    double mSourceCoefficient;

public:

    /**
     * Constructor.
     *
     * @param rCellPopulation reference to the cell population
     * @param duDtCoefficient rate of reaction (defaults to 1.0)
     * @param diffusionCoefficient rate of diffusion (defaults to 1.0)
     * @param sourceCoefficient the source term coefficient (defaults to 0.0)
     */
    CellwiseSourceParabolicPde(AbstractCellPopulation<DIM, DIM>& rCellPopulation,
                               double duDtCoefficient=1.0,
                               double diffusionCoefficient=1.0,
                               double sourceCoefficient=0.0);

    /**
     * @return const reference to the cell population (used in archiving).
     */
    const AbstractCellPopulation<DIM>& rGetCellPopulation() const;

    /**
     * Overridden ComputeDuDtCoefficientFunction() method.
     *
     * @return c
     *
     * @param rX a point in space (unused)
     * @param pdeIndex the index of the PDE (unused)
     */
    virtual double ComputeDuDtCoefficientFunction(const ChastePoint<DIM>& rX,
                                                  unsigned pdeIndex);

    /**
     * Overridden ComputeSourceTerm() method.
     *
     * @return computed source term
     *
     * @param rX the point x at which the source term is computed
     * @param rU the dependent variable, u, at the point x as a vector
     * @param pdeIndex the index of the PDE (unused)
     * @param pElement The mesh element that x is contained in (optional)
     */
    virtual double ComputeSourceTerm(const ChastePoint<DIM>& rX,
                                     c_vector<double,1>& rU,
                                     unsigned pdeIndex,
                                     Element<DIM, DIM>* pElement=nullptr);

    /**
     * Overridden ComputeSourceTermAtNode() method.
     *
     * Note that for cell-wise parabolic PDEs used with CellBasedParabolicPdeSystemSolver
     * this method returns the coefficient of the linear component of the source term.
     *
     * @return computed source term at a node.
     *
     * @param rNode the node at which the source term is computed
     * @param rU the dependent variable, u, at the point x as a vector
     * @param pdeIndex the index of the PDE (unused)
     */
    virtual double ComputeSourceTermAtNode(const Node<DIM>& rNode,
                                           c_vector<double,1>& rU,
                                           unsigned pdeIndex);

    /**
     * Overridden ComputeDiffusionTerm() method.
     *
     * @return computed diffusion term D(x) at a point in space
     *
     * @param rX The point x at which the diffusion term D is computed
     * @param pdeIndex the index of the PDE (unused)
     * @param pElement The mesh element that x is contained in (optional; unused)
     */
    virtual c_matrix<double, DIM, DIM> ComputeDiffusionTerm(const ChastePoint<DIM>& rX,
                                                            unsigned pdeIndex,
                                                            Element<DIM,DIM>* pElement=nullptr);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellwiseSourceParabolicPde)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CellwiseSourceParabolicPde.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const CellwiseSourceParabolicPde<DIM>* t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM, DIM>* p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise a CellwiseSourceParabolicPde.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, CellwiseSourceParabolicPde<DIM>* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM, DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)CellwiseSourceParabolicPde<DIM>(*p_cell_population);
}
}
} // namespace ...

#endif /*CELLWISESOURCEPARABOLICPDE_HPP_*/
