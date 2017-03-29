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

#include "GeneralMonolayerVertexMeshForce.hpp"
#include "MonolayerVertexMeshCustomFunctions.hpp"
#include "VertexBasedCellPopulation.hpp"

GeneralMonolayerVertexMeshForce::GeneralMonolayerVertexMeshForce()
        : AbstractForce<3>(),
          mTargetApicalArea(0),
          mApicalAreaParameter(0),
          mApicalEdgeParameter(0),
          mTargetBasalArea(0),
          mBasalAreaParameter(0),
          mBasalEdgeParameter(0),
          mLateralAreaParameter(0),
          mLateralEdgeParameter(0),
          mTargetVolume(0),
          mVolumeParameter(0)
{
}

GeneralMonolayerVertexMeshForce::~GeneralMonolayerVertexMeshForce()
{
}

c_vector<double, 3> CalculateEdgeGradient(const Node<3>* pNode1, const Node<3>* pNode2)
{
    const c_vector<double, 3>& loc1 = pNode1->rGetLocation();
    const c_vector<double, 3>& loc2 = pNode2->rGetLocation();
    return (loc1 - loc2) / norm_2(loc1 - loc2);
}

void GeneralMonolayerVertexMeshForce::AddForceContribution(AbstractCellPopulation<3>& rCellPopulation)
{
    if (dynamic_cast<VertexBasedCellPopulation<3>*>(&rCellPopulation) == NULL)
    {
        EXCEPTION("GeneralMonolayerVertexMeshForce is to be used with a VertexBasedCellPopulation only"); //LCOV_EXCL_LINE
    }

    // Define some helper variables
    VertexBasedCellPopulation<3>* p_cell_population = static_cast<VertexBasedCellPopulation<3>*>(&rCellPopulation);

    // Add volume contribution
    if (fabs(mVolumeParameter) > DBL_EPSILON)
    {
        this->AddVolumeContribution(p_cell_population);
    }

    if (fabs(mApicalAreaParameter) > DBL_EPSILON || fabs(mBasalAreaParameter) > DBL_EPSILON || fabs(mLateralAreaParameter) > DBL_EPSILON)
    {
        this->AddAreaContribution(p_cell_population);
    }

    if (fabs(mApicalEdgeParameter) > DBL_EPSILON || fabs(mBasalEdgeParameter) > DBL_EPSILON || fabs(mLateralEdgeParameter) > DBL_EPSILON)
    {
        this->AddEdgeContribution(p_cell_population);
    }
}

void GeneralMonolayerVertexMeshForce::AddVolumeContribution(VertexBasedCellPopulation<3>* pCellPopulation)
{
    MutableVertexMesh<3, 3>& rMesh = pCellPopulation->rGetMesh();

    for (unsigned elem_index = 0; elem_index < rMesh.GetNumElements(); ++elem_index)
    {
        const VertexElement<3, 3>* p_elem = rMesh.GetElement(elem_index);
        const double elem_volume = rMesh.GetVolumeOfElement(elem_index);

        for (unsigned local_node_index = 0; local_node_index < p_elem->GetNumNodes(); ++local_node_index)
        {
            Node<3>* p_node = p_elem->GetNode(local_node_index);

            c_vector<double, 3> tmp_v = rMesh.GetVolumeGradientofElementAtNode(p_elem, p_node->GetIndex());
            tmp_v *= -1 * mVolumeParameter * (elem_volume - mTargetVolume);

            p_node->AddAppliedForceContribution(tmp_v);
        }
    }
}

void GeneralMonolayerVertexMeshForce::AddAreaContribution(VertexBasedCellPopulation<3>* pCellPopulation)
{
    MutableVertexMesh<3, 3>& rMesh = pCellPopulation->rGetMesh();

    // Do face contributions.
    for (unsigned face_id = 0; face_id < rMesh.GetNumFaces(); ++face_id)
    {
        const VertexElement<2, 3>* p_face = rMesh.GetFace(face_id);

        switch (GetFaceType(p_face))
        {
            case Monolayer::ApicalValue:
            {
                // Calculate apical area contribution
                if (fabs(mApicalAreaParameter) > 1e-5)
                {
                    const double apical_area = rMesh.CalculateAreaOfFace(p_face);

                    for (unsigned node_id = 0; node_id < p_face->GetNumNodes(); ++node_id)
                    {
                        c_vector<double, 3> result = rMesh.GetAreaGradientOfFaceAtNode(p_face, node_id);
                        result *= -1 * mApicalAreaParameter; // *(apical_area - mTargetApicalArea);
                        p_face->GetNode(node_id)->AddAppliedForceContribution(result);
                    }
                }
                break;
            }
            case Monolayer::BasalValue:
            {
                // Calculate basal area contribution
                if (fabs(mBasalAreaParameter) > 1e-5)
                {
                    const double basal_area = rMesh.CalculateAreaOfFace(p_face);

                    for (unsigned node_id = 0; node_id < p_face->GetNumNodes(); ++node_id)
                    {
                        c_vector<double, 3> result = rMesh.GetAreaGradientOfFaceAtNode(p_face, node_id);
                        result *= -1 * mBasalAreaParameter; // *(basal_area - mTargetBasalArea);
                        p_face->GetNode(node_id)->AddAppliedForceContribution(result);
                    }
                }
                break;
            }
            case Monolayer::LateralValue:
            {
                if (fabs(mLateralAreaParameter) > 1e-5)
                {
                    for (unsigned node_id = 0; node_id < p_face->GetNumNodes(); ++node_id)
                    {
                        c_vector<double, 3> result = rMesh.GetAreaGradientOfFaceAtNode(p_face, node_id);
                        result *= -1 * mLateralAreaParameter;
                        p_face->GetNode(node_id)->AddAppliedForceContribution(result);
                    }
                }
                break;
            }
            default:
                NEVER_REACHED;
        }
    }
}

void GeneralMonolayerVertexMeshForce::AddEdgeContribution(VertexBasedCellPopulation<3>* pCellPopulation)
{
    MutableVertexMesh<3, 3>& rMesh = pCellPopulation->rGetMesh();

    // Do edge contributions.
    for (unsigned face_id = 0; face_id < rMesh.GetNumFaces(); ++face_id)
    {
        const VertexElement<2, 3>* p_face = rMesh.GetFace(face_id);

        if (!IsLateralFace(p_face))
        {
            continue;
        }

        Node<3>* p_node1 = p_face->GetNode(p_face->GetNumNodes() - 1);
        // Calculate apical and basal edge contribution here so that it will be counted once.
        for (unsigned i = 0; i < p_face->GetNumNodes(); ++i)
        {
            Node<3>* p_node2 = p_face->GetNode(i);
            c_vector<double, 3> result = CalculateEdgeGradient(p_node1, p_node2);

            if (GetNodeType(p_node1) == GetNodeType(p_node2))
            {
                if (IsApicalNode(p_node1))
                {
                    result *= -1 * mApicalEdgeParameter;
                }
                else
                {
                    assert(IsBasalNode(p_node1));
                    result *= -1 * mBasalEdgeParameter;
                }
            }
            else
            {
                result *= -1 * mLateralEdgeParameter;
                const std::set<unsigned> s_tmp = GetSharedFaceIndices(p_node1, p_node2);
                result /= s_tmp.size();
            }

            p_node1->AddAppliedForceContribution(result);
            p_node2->AddAppliedForceContribution(-result);

            p_node1 = p_node2;
        }
    }
}

void GeneralMonolayerVertexMeshForce::SetApicalParameters(const double lineParameter, const double areaParameter, const double targetArea)
{
    mTargetApicalArea = targetArea;
    mApicalAreaParameter = areaParameter;
    mApicalEdgeParameter = lineParameter;
}

void GeneralMonolayerVertexMeshForce::SetBasalParameters(const double lineParameter, const double areaParameter, const double targetArea)
{
    mTargetBasalArea = targetArea;
    mBasalAreaParameter = areaParameter;
    mBasalEdgeParameter = lineParameter;
}

void GeneralMonolayerVertexMeshForce::SetLateralParameter(const double lineParameter, const double areaParameter)
{
    mLateralEdgeParameter = lineParameter;
    mLateralAreaParameter = areaParameter;
}

void GeneralMonolayerVertexMeshForce::SetVolumeParameters(const double volumeParameter, const double targetVolume)
{
    mTargetVolume = targetVolume;
    mVolumeParameter = volumeParameter;
}

void GeneralMonolayerVertexMeshForce::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<TargetApicalArea>" << mTargetApicalArea << "</TargetApicalArea>\n";
    *rParamsFile << "\t\t\t<ApicalAreaParameter>" << mApicalAreaParameter << "</ApicalAreaParameter>\n";
    *rParamsFile << "\t\t\t<ApicalEdgeParameter>" << mApicalEdgeParameter << "</ApicalEdgeParameter>\n";

    *rParamsFile << "\t\t\t<TargetBasalArea>" << mTargetBasalArea << "</TargetBasalArea>\n";
    *rParamsFile << "\t\t\t<BasalAreaParameter>" << mBasalAreaParameter << "</BasalAreaParameter>\n";
    *rParamsFile << "\t\t\t<BasalEdgeParameter>" << mBasalEdgeParameter << "</BasalEdgeParameter>\n";
    *rParamsFile << "\t\t\t<LateralAreaParameter>" << mLateralAreaParameter << "</LateralAreaParameter>\n";
    *rParamsFile << "\t\t\t<LateralEdgeParameter>" << mLateralEdgeParameter << "</LateralEdgeParameter>\n";
    *rParamsFile << "\t\t\t<TargetVolume>" << mTargetVolume << "</TargetVolume>\n";
    *rParamsFile << "\t\t\t<VolumeParameter>" << mVolumeParameter << "</VolumeParameter>\n";

    AbstractForce<3>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(GeneralMonolayerVertexMeshForce)
