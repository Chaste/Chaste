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

#include "BidomainTissue.hpp"

#include "DistributedVector.hpp"
#include "AxisymmetricConductivityTensors.hpp"
#include "OrthotropicConductivityTensors.hpp"
#include "ChastePoint.hpp"
#include "AbstractChasteRegion.hpp"

template <unsigned SPACE_DIM>
BidomainTissue<SPACE_DIM>::BidomainTissue(
            AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory,
            bool exchangeHalos)
    : AbstractCardiacTissue<SPACE_DIM>(pCellFactory, exchangeHalos)
{
    CreateExtracellularConductivityTensors();
}

template <unsigned SPACE_DIM>
BidomainTissue<SPACE_DIM>::BidomainTissue(AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>* pMesh)
        :  AbstractCardiacTissue<SPACE_DIM>(pMesh)
{
    CreateExtracellularConductivityTensors();
}

template <unsigned SPACE_DIM>
void BidomainTissue<SPACE_DIM>::CreateExtracellularConductivityTensors()
{
    if (this->mpConfig->IsMeshProvided() && this->mpConfig->GetLoadMesh())
    {
        assert(this->mFibreFilePathNoExtension != "");

        switch (this->mpConfig->GetConductivityMedia())
        {
            case cp::media_type::Orthotropic:
            {
                mpExtracellularConductivityTensors =  new OrthotropicConductivityTensors<SPACE_DIM,SPACE_DIM>;
                FileFinder ortho_file(this->mFibreFilePathNoExtension + ".ortho", RelativeTo::AbsoluteOrCwd);
                assert(ortho_file.Exists());
                mpExtracellularConductivityTensors->SetFibreOrientationFile(ortho_file);
                break;
            }

            case cp::media_type::Axisymmetric:
            {
                mpExtracellularConductivityTensors =  new AxisymmetricConductivityTensors<SPACE_DIM,SPACE_DIM>;
                FileFinder axi_file(this->mFibreFilePathNoExtension + ".axi", RelativeTo::AbsoluteOrCwd);
                assert(axi_file.Exists());
                mpExtracellularConductivityTensors->SetFibreOrientationFile(axi_file);
                break;
            }

            case cp::media_type::NoFibreOrientation:
                mpExtracellularConductivityTensors =  new OrthotropicConductivityTensors<SPACE_DIM,SPACE_DIM>;
                break;

            default :
                NEVER_REACHED;
        }
    }
    else // Slab defined in config file or SetMesh() called; no fibre orientation assumed
    {
        mpExtracellularConductivityTensors =  new OrthotropicConductivityTensors<SPACE_DIM,SPACE_DIM>;
    }

    c_vector<double, SPACE_DIM> extra_conductivities;
    this->mpConfig->GetExtracellularConductivities(extra_conductivities);

    // this definition must be here (and not inside the if statement) because SetNonConstantConductivities() will keep
    // a pointer to it and we don't want it to go out of scope before Init() is called
    unsigned num_local_elements = this->mpMesh->GetNumLocalElements();
    std::vector<c_vector<double, SPACE_DIM> > hetero_extra_conductivities;

    if (this->mpConfig->GetConductivityHeterogeneitiesProvided())
    {
        try
        {
            assert(hetero_extra_conductivities.size()==0);
            //initialise with the values of the default conductivity tensor
            hetero_extra_conductivities.resize(num_local_elements, extra_conductivities);
        }
        // LCOV_EXCL_START
        catch(std::bad_alloc &r_bad_alloc)
        {
            std::cout << "Failed to allocate std::vector of size " << num_local_elements << std::endl;
            PetscTools::ReplicateException(true);
            throw r_bad_alloc;
        }
        // LCOV_EXCL_STOP

        PetscTools::ReplicateException(false);

        std::vector<boost::shared_ptr<AbstractChasteRegion<SPACE_DIM> > > conductivities_heterogeneity_areas;
        std::vector< c_vector<double,3> > intra_h_conductivities;
        std::vector< c_vector<double,3> > extra_h_conductivities;
        HeartConfig::Instance()->GetConductivityHeterogeneities(conductivities_heterogeneity_areas,
                                                                intra_h_conductivities,
                                                                extra_h_conductivities);

        unsigned local_element_index = 0;

        for (typename AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>::ElementIterator iter = (this->mpMesh)->GetElementIteratorBegin();
             iter != (this->mpMesh)->GetElementIteratorEnd();
             ++iter)
        {
            // If element centroid is contained in the region
            ChastePoint<SPACE_DIM> element_centroid(iter->CalculateCentroid());
            for (unsigned region_index=0; region_index< conductivities_heterogeneity_areas.size(); region_index++)
            {
                // If element centroid is contained in the region
                if (conductivities_heterogeneity_areas[region_index]->DoesContain(element_centroid))
                {
                    //We don't use ublas vector assignment here, because we might be getting a subvector of a 3-vector
                    for (unsigned i=0; i<SPACE_DIM; i++)
                    {
                        hetero_extra_conductivities[local_element_index][i] = extra_h_conductivities[region_index][i];
                    }
                }
            }
            local_element_index++;
        }
        mpExtracellularConductivityTensors->SetNonConstantConductivities(&hetero_extra_conductivities);
    }
    else
    {
        mpExtracellularConductivityTensors->SetConstantConductivities(extra_conductivities);
    }

    mpExtracellularConductivityTensors->Init(this->mpMesh);
}

template <unsigned SPACE_DIM>
BidomainTissue<SPACE_DIM>::~BidomainTissue()
{
    if (mpExtracellularConductivityTensors)
    {
        delete mpExtracellularConductivityTensors;
    }
}


template <unsigned SPACE_DIM>
const c_matrix<double, SPACE_DIM, SPACE_DIM>& BidomainTissue<SPACE_DIM>::rGetExtracellularConductivityTensor(unsigned elementIndex)
{
    assert(mpExtracellularConductivityTensors);
    if (this->mpConductivityModifier==NULL)
    {
        return (*mpExtracellularConductivityTensors)[elementIndex];
    }
    else
    {
        return this->mpConductivityModifier->rGetModifiedConductivityTensor(elementIndex, (*mpExtracellularConductivityTensors)[elementIndex], 1u);
    }
}

// Explicit instantiation
template class BidomainTissue<1>;
template class BidomainTissue<2>;
template class BidomainTissue<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BidomainTissue)
