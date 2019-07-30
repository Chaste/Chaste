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

#include "MonodomainProblem.hpp"

#include "Exception.hpp"
#include "ReplicatableVector.hpp"
#include "MonodomainSolver.hpp"
#include "OperatorSplittingMonodomainSolver.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>* MonodomainProblem<ELEMENT_DIM, SPACE_DIM>::CreateCardiacTissue()
{
    mpMonodomainTissue = new MonodomainTissue<ELEMENT_DIM,SPACE_DIM>(this->mpCellFactory, HeartConfig::Instance()->GetUseStateVariableInterpolation());
    return mpMonodomainTissue;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, 1>* MonodomainProblem<ELEMENT_DIM, SPACE_DIM>::CreateSolver()
{
    assert(mpMonodomainTissue);
    /*
     * NOTE: The this->mpBoundaryConditionsContainer.get() lines below convert a
     * boost::shared_ptr to a normal pointer, as this is what the solvers are
     * expecting. We have to be a bit careful though as boost could decide to delete
     * them whenever it feels like as it won't count the assembers as using them.
     *
     * As long as they are kept as member variables here for as long as they are
     * required in the solvers it should all work OK.
     */

    if (HeartConfig::Instance()->GetUseReactionDiffusionOperatorSplitting())
    {
        return new OperatorSplittingMonodomainSolver<ELEMENT_DIM,SPACE_DIM>(this->mpMesh,
                                                                            mpMonodomainTissue,
                                                                            this->mpBoundaryConditionsContainer.get());
    }
    else
    {
        return new MonodomainSolver<ELEMENT_DIM,SPACE_DIM>(this->mpMesh,
                                                           mpMonodomainTissue,
                                                           this->mpBoundaryConditionsContainer.get());
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainProblem<ELEMENT_DIM, SPACE_DIM>::MonodomainProblem(AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>* pCellFactory)
        : AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, 1>(pCellFactory),
          mpMonodomainTissue(NULL)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainProblem<ELEMENT_DIM, SPACE_DIM>::MonodomainProblem()
    : AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, 1>(),
      mpMonodomainTissue(NULL)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainProblem<ELEMENT_DIM, SPACE_DIM>::~MonodomainProblem()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainTissue<ELEMENT_DIM,SPACE_DIM> * MonodomainProblem<ELEMENT_DIM, SPACE_DIM>::GetMonodomainTissue()
{
    assert(mpMonodomainTissue != NULL);
    return mpMonodomainTissue;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainProblem<ELEMENT_DIM, SPACE_DIM>::WriteInfo(double time)
{
    if (PetscTools::AmMaster())
    {
        std::cout << "Solved to time " << time << "\n" << std::flush;
    }

    double v_max, v_min;

    VecMax( this->mSolution, PETSC_NULL, &v_max );
    VecMin( this->mSolution, PETSC_NULL, &v_min );

    if (PetscTools::AmMaster())
    {
        std::cout << " V = " << "[" <<v_min << ", " << v_max << "]" << "\n" << std::flush;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainProblem<ELEMENT_DIM, SPACE_DIM>::DefineWriterColumns(bool extending)
{
    AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,1>::DefineWriterColumns(extending);
    AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,1>::DefineExtraVariablesWriterColumns(extending);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainProblem<ELEMENT_DIM, SPACE_DIM>::WriteOneStep(double time, Vec voltageVec)
{
    this->mpWriter->PutUnlimitedVariable(time);
    this->mpWriter->PutVector(this->mVoltageColumnId, voltageVec);
    AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,1>::WriteExtraVariablesOneStep();
}

// Explicit instantiation
template class MonodomainProblem<1,1>;
template class MonodomainProblem<1,2>;
template class MonodomainProblem<1,3>;
template class MonodomainProblem<2,2>;
template class MonodomainProblem<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS2(MonodomainProblem, 1, 1)
EXPORT_TEMPLATE_CLASS2(MonodomainProblem, 1, 2)
EXPORT_TEMPLATE_CLASS2(MonodomainProblem, 1, 3)
EXPORT_TEMPLATE_CLASS2(MonodomainProblem, 2, 2)
EXPORT_TEMPLATE_CLASS2(MonodomainProblem, 3, 3)
