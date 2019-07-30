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

#include "MonodomainPurkinjeProblem.hpp"
#include "MonodomainPurkinjeSolver.hpp"
#include "Exception.hpp"
#include "ReplicatableVector.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>* MonodomainPurkinjeProblem<ELEMENT_DIM, SPACE_DIM>::CreateCardiacTissue()
{
    return new MonodomainTissue<ELEMENT_DIM,SPACE_DIM>(this->mpCellFactory, HeartConfig::Instance()->GetUseStateVariableInterpolation());
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, 2>* MonodomainPurkinjeProblem<ELEMENT_DIM, SPACE_DIM>::CreateSolver()
{
    assert(!HeartConfig::Instance()->GetUseReactionDiffusionOperatorSplitting());
    MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>* p_mesh = dynamic_cast<MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>*>(this->mpMesh);
    EXCEPT_IF_NOT(p_mesh); ///\todo #2017 give a nice error
    MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* p_tissue = dynamic_cast<MonodomainTissue<ELEMENT_DIM,SPACE_DIM>*>(this->GetTissue());
    assert(p_tissue);
    return new MonodomainPurkinjeSolver<ELEMENT_DIM,SPACE_DIM>(p_mesh,
                                                               p_tissue,
                                                               this->mpBoundaryConditionsContainer.get());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainPurkinjeProblem<ELEMENT_DIM,SPACE_DIM>::CreateMeshFromHeartConfig()
{
    this->mpMesh = new MixedDimensionMesh<ELEMENT_DIM, SPACE_DIM>(HeartConfig::Instance()->GetMeshPartitioning());
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Vec MonodomainPurkinjeProblem<ELEMENT_DIM, SPACE_DIM>::CreateInitialCondition()
{
    Vec init_cond = AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,2>::CreateInitialCondition();

    // Get the Purkinje voltage stripe
    DistributedVector ic = this->mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(init_cond);
    DistributedVector::Stripe purkinje_voltage_stripe = DistributedVector::Stripe(ic, 1);

    for (DistributedVector::Iterator index = ic.Begin();
         index != ic.End();
         ++index)
    {
        purkinje_voltage_stripe[index] = this->GetTissue()->GetPurkinjeCell(index.Global)->GetVoltage();
        // Note that it doesn't matter if the cell is fake, as this will give a zero voltage for nodes not in the
        // Purkinje network.
    }
    ic.Restore();

    return init_cond;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainPurkinjeProblem<ELEMENT_DIM, SPACE_DIM>::MonodomainPurkinjeProblem(AbstractPurkinjeCellFactory<ELEMENT_DIM,SPACE_DIM>* pCellFactory)
        : AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, 2>(pCellFactory),
          mPurkinjeVoltageColumnId(UNSIGNED_UNSET)
{
}


// LCOV_EXCL_START
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainPurkinjeProblem<ELEMENT_DIM, SPACE_DIM>::MonodomainPurkinjeProblem()
    : AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, 2>(),
      mPurkinjeVoltageColumnId(UNSIGNED_UNSET)
{
}
// LCOV_EXCL_STOP

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonodomainPurkinjeProblem<ELEMENT_DIM, SPACE_DIM>::~MonodomainPurkinjeProblem()
{
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainPurkinjeProblem<ELEMENT_DIM, SPACE_DIM>::WriteInfo(double time)
{
    if (PetscTools::AmMaster())
    {
        std::cout << "Solved to time " << time << "\n" << std::flush;
    }

    double v_max, v_min, v_purk_max, v_purk_min;

    VecStrideMax( this->mSolution, 0, PETSC_NULL, &v_max );
    VecStrideMin( this->mSolution, 0, PETSC_NULL, &v_min );

    VecStrideMax( this->mSolution, 1, PETSC_NULL, &v_purk_max );
    VecStrideMin( this->mSolution, 1, PETSC_NULL, &v_purk_min );

    // avoid printing 1e-320 etc
    if(fabs(v_purk_min)<1e-20)
    {
        // It's very hard to hit this, as you'd need the whole Purkinje system to be activated
        v_purk_min = 0.0; // LCOV_EXCL_LINE
    }

    if(fabs(v_purk_max)<1e-20)
    {
        v_purk_max = 0.0;
    }

    if (PetscTools::AmMaster())
    {
        std::cout << " V_myo; V_purk = " << "[" <<v_min << ", " << v_max << "]" << ";\t"
                  << "[" << v_purk_min << ", " << v_purk_max << "]" << "\n"
                  << std::flush;
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainPurkinjeProblem<ELEMENT_DIM, SPACE_DIM>::DefineWriterColumns(bool extending)
{
    AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM,2>::DefineWriterColumns(extending);
    if (extending)
    {
        mPurkinjeVoltageColumnId = this->mpWriter->GetVariableByName("V_purk");
    }
    else
    {
        mPurkinjeVoltageColumnId = this->mpWriter->DefineVariable("V_purk","mV");
    }
    AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM,2>::DefineExtraVariablesWriterColumns(extending);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonodomainPurkinjeProblem<ELEMENT_DIM, SPACE_DIM>::WriteOneStep(double time, Vec voltageVec)
{
    this->mpWriter->PutUnlimitedVariable(time);
    std::vector<int> variable_ids;
    variable_ids.push_back(this->mVoltageColumnId);
    variable_ids.push_back(mPurkinjeVoltageColumnId);
    this->mpWriter->PutStripedVector(variable_ids, voltageVec);
    AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,2>::WriteExtraVariablesOneStep();
}


/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class MonodomainPurkinjeProblem<2,2>;
template class MonodomainPurkinjeProblem<3,3>;



// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS2(MonodomainPurkinjeProblem, 2, 2)
EXPORT_TEMPLATE_CLASS2(MonodomainPurkinjeProblem, 3, 3)
