/*

Copyright (C) University of Oxford, 2005-2012

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
        catch(std::bad_alloc &r_bad_alloc)
        {
#define COVERAGE_IGNORE
            std::cout << "Failed to allocate std::vector of size " << num_local_elements << std::endl;
            PetscTools::ReplicateException(true);
            throw r_bad_alloc;
#undef COVERAGE_IGNORE
        }
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
            // if element centroid is contained in the region
            ChastePoint<SPACE_DIM> element_centroid(iter->CalculateCentroid());
            for (unsigned region_index=0; region_index< conductivities_heterogeneity_areas.size(); region_index++)
            {
                // if element centroid is contained in the region
              if ( conductivities_heterogeneity_areas[region_index]->DoesContain( element_centroid ) )
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
    if(this->mpConductivityModifier==NULL)
    {
        return (*mpExtracellularConductivityTensors)[elementIndex];
    }
    else
    {
        return this->mpConductivityModifier->rGetModifiedConductivityTensor(elementIndex, (*mpExtracellularConductivityTensors)[elementIndex]);
    }
}


/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class BidomainTissue<1>;
template class BidomainTissue<2>;
template class BidomainTissue<3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BidomainTissue)
