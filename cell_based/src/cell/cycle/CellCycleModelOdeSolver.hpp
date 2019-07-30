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

#ifndef CELLCYCLEMODELODESOLVER_HPP_
#define CELLCYCLEMODELODESOLVER_HPP_

#include <boost/utility.hpp>

#include "ChasteSerialization.hpp"

#include "AbstractCellCycleModelOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "CvodeAdaptor.hpp"

/**
 * A concrete implementation of AbstractCellCycleModelOdeSolver, that uses templates
 * to provide an implementation for any pair of cell-cycle model and ODE solver classes.
 *
 * All ODE-based cell-cycle model developers need to do is set mpOdeSolver in their constructor:
 *   mpOdeSolver = CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER>::Instance();
 *
 * This class contains all the machinery to make it a singleton, hence providing
 * exactly one instance per pair of values of the template parameters.
 */
template <class CELL_CYCLE_MODEL, class ODE_SOLVER>
class CellCycleModelOdeSolver : public AbstractCellCycleModelOdeSolver, private boost::noncopyable
{
private:
    /** The single instance of this class, for this ODE_SOLVER. */
    static boost::shared_ptr<CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER> > mpInstance;

    /** Default constructor. Not user accessible; to obtain an instance of this class use the Instance method. */
    CellCycleModelOdeSolver();

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModelOdeSolver>(*this);
        archive & mpInstance;
    }

public:
    /** @return a pointer to the singleton instance, creating it if necessary. */
    static boost::shared_ptr<CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER> > Instance();

    /** @return whether the instance in existence and fully set up. */
    bool IsSetUp();

    /** Initialise the ODE solver. */
    void Initialise();

    /**
     * @return true iff this is an adaptive solver such as CVODE for which it is safe to set the 'timestep'
     * to be the outer simulation timestep, because the ODE solver will use this as its maximum, not actual,
     * timestep.
     *
     * By default calls the base class version; it is defined here so that specializations can override it.
     */
    virtual bool IsAdaptive();

    /** Reset the instance. */
    void Reset();
};

/** Definition of the instance static member. */
template<class CELL_CYCLE_MODEL, class ODE_SOLVER>
boost::shared_ptr<CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER> > CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER>::mpInstance;


template<class CELL_CYCLE_MODEL, class ODE_SOLVER>
CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER>::CellCycleModelOdeSolver()
    : AbstractCellCycleModelOdeSolver()
{
    /**
     * The semantics of shared_ptr are different from normal pointers; we don't
     * care if a second instance is constructed when loading an archive, since
     * archiving the shared_ptr will ensure that the 'singleton' is correctly
     * serialized. Thus, here we do not require an assertion that mpInstance is
     * NULL, as we do in the constructors of the SimulationTime singleton class.
     */
}

template<class CELL_CYCLE_MODEL, class ODE_SOLVER>
boost::shared_ptr<CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER> > CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER>::Instance()
{
    if (!mpInstance)
    {
        mpInstance.reset(new CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER>);
    }
    return mpInstance;
}

template<class CELL_CYCLE_MODEL, class ODE_SOLVER>
bool CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER>::IsSetUp()
{
    return static_cast<bool>(mpOdeSolver.get());
}

template<class CELL_CYCLE_MODEL, class ODE_SOLVER>
void CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER>::Initialise()
{
    mpOdeSolver.reset(new ODE_SOLVER);
    // If this is a CVODE solver we need to tell it to reset. Otherwise
    // the fact this is a singleton will lead to all sorts of problems
    // as CVODE will have the internal state for the wrong ODE system!
#ifdef CHASTE_CVODE
    if (boost::dynamic_pointer_cast<CvodeAdaptor>(mpOdeSolver))
    {
        (boost::static_pointer_cast<CvodeAdaptor>(mpOdeSolver))->SetForceReset(true);
    }
#endif //CHASTE_CVODE
}

template<class CELL_CYCLE_MODEL, class ODE_SOLVER>
bool CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER>::IsAdaptive()
{
    return AbstractCellCycleModelOdeSolver::IsAdaptive();
}

template<class CELL_CYCLE_MODEL, class ODE_SOLVER>
void CellCycleModelOdeSolver<CELL_CYCLE_MODEL, ODE_SOLVER>::Reset()
{
    ///\todo Consider whether Reset() could be moved to the abstract class
}

/**
 * Specialization for BackwardEulerIvpOdeSolver, whose constructor requires
 * an argument.
 */
template<class CELL_CYCLE_MODEL>
class CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver> : public AbstractCellCycleModelOdeSolver, private boost::noncopyable
{
private:
    /** The single instance of this class, for this ODE_SOLVER. */
    static boost::shared_ptr<CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver> > mpInstance;

    /** Default constructor. Not user accessible; to obtain an instance of this class use the Instance method. */
    CellCycleModelOdeSolver();

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModelOdeSolver>(*this);
        archive & mpInstance;
    }

public:
    /** @return a pointer to the singleton instance, creating it if necessary. */
    static boost::shared_ptr<CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver> > Instance();

    /** @return whether the instance in existence and fully set up. */
    bool IsSetUp();

    /** Initialise the ODE solver. */
    void Initialise();

    /** Reset the instance. */
    void Reset();
};

template<class CELL_CYCLE_MODEL>
boost::shared_ptr<CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver> > CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver>::mpInstance;

template<class CELL_CYCLE_MODEL>
CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver>::CellCycleModelOdeSolver()
    : AbstractCellCycleModelOdeSolver()
{
}

template<class CELL_CYCLE_MODEL>
boost::shared_ptr<CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver> > CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver>::Instance()
{
    if (!mpInstance)
    {
        mpInstance.reset(new CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver>);
    }
    return mpInstance;
}

template<class CELL_CYCLE_MODEL>
bool CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver>::IsSetUp()
{
    return mpOdeSolver && (mSizeOfOdeSystem != UNSIGNED_UNSET);
}

template<class CELL_CYCLE_MODEL>
void CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver>::Initialise()
{
    if (mSizeOfOdeSystem == UNSIGNED_UNSET)
    {
        EXCEPTION("SetSizeOfOdeSystem() must be called before calling Initialise()");
    }
    mpOdeSolver.reset(new BackwardEulerIvpOdeSolver(mSizeOfOdeSystem));
}

template<class CELL_CYCLE_MODEL>
void CellCycleModelOdeSolver<CELL_CYCLE_MODEL, BackwardEulerIvpOdeSolver>::Reset()
{
    mSizeOfOdeSystem = UNSIGNED_UNSET;
    mpOdeSolver.reset();
}

#endif /*CELLCYCLEMODELODESOLVER_HPP_*/
