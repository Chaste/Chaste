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


#include "AbstractPurkinjeCellFactory.hpp"
#include "PurkinjeVentricularJunctionStimulus.hpp"
#include "MultiStimulus.hpp"
#include "HeartConfig.hpp"
#include "Warnings.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractPurkinjeCellFactory<ELEMENT_DIM,SPACE_DIM>::AbstractPurkinjeCellFactory()
    : AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>(),
      mpMixedDimensionMesh(NULL)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractPurkinjeCellFactory<ELEMENT_DIM,SPACE_DIM>::ReadJunctionsFile()
{
    std::string pvj_file_name;
    bool file_specified = true;

    try
    {
        // HeartConfig::Instance()->GetMeshName() will throw an exception if no mesh name is defined
        pvj_file_name = HeartConfig::Instance()->GetMeshName() + ".pvj";
    }
    catch(Exception&)
    {
        file_specified = false;
    }

    FileFinder junction_file(pvj_file_name, RelativeTo::AbsoluteOrCwd);
    if (!file_specified || !junction_file.Exists() )
    {
        // In this case we expect the user to specify PVJ programmatically, so return immediately.
        WARNING("No Purkinje-Ventricular junction (.pvj) file found. Junctions must be specified manually.");
        return;
    }

    std::ifstream junction_stream(junction_file.GetAbsolutePath().c_str());

    if (!junction_stream.good())
    {   // file couldn't be opened
        EXCEPTION("Couldn't open data file: " << junction_file.GetAbsolutePath());
    }

    // Reads in file defining nodes and resistance (separated by space)
    while (junction_stream.good())
    {
        std::string this_line;
        getline(junction_stream, this_line);

        if (this_line=="" || this_line=="\r")
        {
            if (junction_stream.eof())
            {   // If the blank line is the last line carry on OK.
                break;
            }
            else
            {
                EXCEPTION("No data found on line in file: " << junction_file.GetAbsolutePath());
            }
        }
        std::stringstream line(this_line);

        unsigned node_id;
        line >> node_id;
        double resistance;
        line >> resistance;

        if (mpMixedDimensionMesh->rGetNodePermutation().size() != 0) //Do we have a permuted mesh?
        {
            unsigned mapped_node_id = mpMixedDimensionMesh->rGetNodePermutation()[node_id];

            mJunctionMap[mapped_node_id] = resistance;
        }
        else
        {
            mJunctionMap[node_id] = resistance;
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractPurkinjeCellFactory<ELEMENT_DIM,SPACE_DIM>::CreateJunction(const Node<SPACE_DIM>* pNode,
                                                                        AbstractCardiacCellInterface* pPurkinjeCell,
                                                                        AbstractCardiacCellInterface* pCardiacCell,
                                                                        double resistance)
{
    // Figure out the effective resistance for this mesh, in kOhm.cm^3
    if (pNode) // Should always be provided apart from low-level tests!
    {
        assert(mpMixedDimensionMesh);
        typedef typename MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>::CableRangeAtNode CableRangeAtNode;
        CableRangeAtNode cable_range = mpMixedDimensionMesh->GetCablesAtNode(pNode);
        double total_cross_sectional_area = 0.0;
        for (typename MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>::NodeCableIterator iter=cable_range.first;
             iter != cable_range.second;
             ++iter)
        {
            Element<1u,SPACE_DIM>* p_cable = iter->second;
            double cable_radius = p_cable->GetAttribute();
            total_cross_sectional_area += M_PI*cable_radius*cable_radius;
        }
        resistance *= total_cross_sectional_area / HeartConfig::Instance()->GetPurkinjeSurfaceAreaToVolumeRatio();
    }
    // Create the junction stimuli, and associate them with the cells
    boost::shared_ptr<PurkinjeVentricularJunctionStimulus> p_pvj_ventricular_stim(new PurkinjeVentricularJunctionStimulus(resistance));
    boost::shared_ptr<PurkinjeVentricularJunctionStimulus> p_pvj_purkinje_stim(new PurkinjeVentricularJunctionStimulus(resistance));
    p_pvj_purkinje_stim->SetAppliedToPurkinjeCellModel();
    p_pvj_ventricular_stim->SetVentricularCellModel(pCardiacCell);
    p_pvj_ventricular_stim->SetPurkinjeCellModel(pPurkinjeCell);
    p_pvj_purkinje_stim->SetVentricularCellModel(pCardiacCell);
    p_pvj_purkinje_stim->SetPurkinjeCellModel(pPurkinjeCell);

    // Create new combined stimuli which add the junction stimuli to those already in the cells
    boost::shared_ptr<MultiStimulus> p_multi_stim_ventricular(new MultiStimulus);
    p_multi_stim_ventricular->AddStimulus(p_pvj_ventricular_stim);
    p_multi_stim_ventricular->AddStimulus(pCardiacCell->GetStimulusFunction());
    pCardiacCell->SetStimulusFunction(p_multi_stim_ventricular);

    boost::shared_ptr<MultiStimulus> p_multi_stim_purkinje(new MultiStimulus);
    p_multi_stim_purkinje->AddStimulus(p_pvj_purkinje_stim);
    p_multi_stim_purkinje->AddStimulus(pPurkinjeCell->GetStimulusFunction());
    pPurkinjeCell->SetStimulusFunction(p_multi_stim_purkinje);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractPurkinjeCellFactory<ELEMENT_DIM,SPACE_DIM>::SetMesh(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
{
    mpMixedDimensionMesh = dynamic_cast<MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>*>(pMesh);
    if (mpMixedDimensionMesh == NULL)
    {
        EXCEPTION("AbstractPurkinjeCellFactory must take a MixedDimensionMesh");
    }
    mLocalPurkinjeNodes.clear();
    for (typename MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>::CableElementIterator iter = mpMixedDimensionMesh->GetCableElementIteratorBegin();
          iter != mpMixedDimensionMesh->GetCableElementIteratorEnd();
          ++iter)
    {
        mLocalPurkinjeNodes.insert((*iter)->GetNodeGlobalIndex(0u));
        mLocalPurkinjeNodes.insert((*iter)->GetNodeGlobalIndex(1u));
    }
    AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>::SetMesh(pMesh);

    ReadJunctionsFile();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCardiacCellInterface*  AbstractPurkinjeCellFactory<ELEMENT_DIM,SPACE_DIM>::CreatePurkinjeCellForNode(
        Node<SPACE_DIM>* pNode,
        AbstractCardiacCellInterface* pCardiacCell)
{
    unsigned node_index = pNode->GetIndex();
    if (mLocalPurkinjeNodes.count(node_index)>0)
    {
        return CreatePurkinjeCellForTissueNode(pNode, pCardiacCell);
    }
    else
    {
        return new FakeBathCell(this->mpSolver, this->mpZeroStimulus);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>* AbstractPurkinjeCellFactory<ELEMENT_DIM,SPACE_DIM>::GetMixedDimensionMesh()
{
    if (mpMixedDimensionMesh == NULL)
    {
        EXCEPTION("The mixed dimension mesh object has not been set in the cell factory");
    }
    return mpMixedDimensionMesh;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractPurkinjeCellFactory<ELEMENT_DIM,SPACE_DIM>::CreateJunctionFromFile(const Node<SPACE_DIM>* pNode,
                                                                                AbstractCardiacCellInterface* pPurkinjeCell,
                                                                                AbstractCardiacCellInterface* pCardiacCell)
{
    std::map<unsigned, double>::iterator iter = mJunctionMap.find(pNode->GetIndex());
    if (iter != mJunctionMap.end())
    {
        CreateJunction(pNode, pPurkinjeCell, pCardiacCell, iter->second);
    }
}

// Explicit instantiation
template class AbstractPurkinjeCellFactory<1,1>;
template class AbstractPurkinjeCellFactory<2,2>;
template class AbstractPurkinjeCellFactory<3,3>;
template class AbstractPurkinjeCellFactory<1,2>;
template class AbstractPurkinjeCellFactory<1,3>;
