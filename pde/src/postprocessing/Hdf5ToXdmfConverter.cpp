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

#include "Hdf5ToXdmfConverter.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Hdf5ToXdmfConverter<ELEMENT_DIM, SPACE_DIM>::Hdf5ToXdmfConverter(const FileFinder& rInputDirectory,
                                                                 const std::string& rFileBaseName,
                                                                 AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh)
    : AbstractHdf5Converter<ELEMENT_DIM,SPACE_DIM>(rInputDirectory, rFileBaseName, pMesh, "xdmf_output", 0u),
      XdmfMeshWriter<ELEMENT_DIM,SPACE_DIM>(rInputDirectory.GetRelativePath(FileFinder("", RelativeTo::ChasteTestOutput)) + "/xdmf_output",
                                            rFileBaseName, false /* Not cleaning directory*/)
{
    // Set number of time steps
    std::vector<double> time_values = this->mpReader->GetUnlimitedDimensionValues();
    unsigned num_timesteps = time_values.size();
    this->mNumberOfTimePoints = num_timesteps;
    // Set time step size
    if (num_timesteps > 1)
    {
        this->mTimeStep = time_values[1] - time_values[0];
    }
    // Write
    this->WriteFilesUsingMesh(*pMesh);
}

#ifndef _MSC_VER

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void Hdf5ToXdmfConverter<ELEMENT_DIM, SPACE_DIM>::AddDataOnNodes(XERCES_CPP_NAMESPACE_QUALIFIER DOMElement* pGridElement,
                                                                 XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument* pDomDocument,
                                                                 unsigned timeStep)
{
    // Use xerces namespace for convenience
    XERCES_CPP_NAMESPACE_USE

    unsigned num_timesteps = this->mpReader->GetUnlimitedDimensionValues().size();

    // Loop over variables
    for (unsigned var_index=0; var_index<this->mNumVariables; var_index++)
    {
        std::string variable_name = this->mpReader->GetVariableNames()[var_index];

        /*
         * e.g. <Attribute Center="Node" Name="V">
         */
        DOMElement* p_attr_element =  pDomDocument->createElement(X("Attribute"));
        p_attr_element->setAttribute(X("Name"), X(variable_name));
        p_attr_element->setAttribute(X("Center"), X("Node"));
        pGridElement->appendChild(p_attr_element);

        /*
         * e.g. <DataItem Dimensions="1 12 1" ItemType="HyperSlab">
         */
        DOMElement* p_hype_element =  pDomDocument->createElement(X("DataItem"));
        p_hype_element->setAttribute(X("ItemType"), X("HyperSlab"));
        std::stringstream dim_stream;

        // First index is time value, second is number of nodes, third is variable index
        ///\todo #1157 Make this work in parallel
        /* DistributedVectorFactory* p_factory = AbstractHdf5Converter<ELEMENT_DIM, SPACE_DIM>::mpMesh->GetDistributedVectorFactory();
        dim_stream << "1 " << p_factory->GetHigh()-p_factory->GetLow() << " 1"; */
        unsigned num_nodes = AbstractHdf5Converter<ELEMENT_DIM, SPACE_DIM>::mpMesh->GetNumNodes();
        dim_stream << "1 " << num_nodes << " 1";
        p_hype_element->setAttribute(X("Dimensions"), X(dim_stream.str()));
        p_attr_element->appendChild(p_hype_element);

        /*
         * e.g. <DataItem Dimensions="3 3" Format="XML">
         */
        DOMElement* p_xml_element =  pDomDocument->createElement(X("DataItem"));
        p_xml_element->setAttribute(X("Format"), X("XML"));
        p_xml_element->setAttribute(X("Dimensions"), X("3 3"));
        p_hype_element->appendChild(p_xml_element);

        /*
         * e.g. 0 0 0 1 1 1 1 12 1
         */
        std::stringstream XMLStream;
        XMLStream << timeStep << " 0 " << var_index << " ";
        XMLStream << "1 1 1 ";
        /* XMLStream << "1 " << p_factory->GetHigh()-p_factory->GetLow() << " 1"; */
        XMLStream << "1 " << num_nodes << " 1";
        DOMText* p_xml_text = pDomDocument->createTextNode(X(XMLStream.str()));
        p_xml_element->appendChild(p_xml_text);

        /*
         * e.g. <DataItem Dimensions="2 12 2" Format="HDF" NumberType="Float" Precision="8">
         */
        DOMElement* p_hdf_element =  pDomDocument->createElement(X("DataItem"));
        p_hdf_element->setAttribute(X("Format"), X("HDF"));
        p_hdf_element->setAttribute(X("NumberType"), X("Float"));
        p_hdf_element->setAttribute(X("Precision"), X("8"));
        std::stringstream hdf_dims_stream;
        /* hdf_dims_stream << num_timesteps << " " << p_factory->GetHigh()-p_factory->GetLow() << " " << this->mNumVariables; */
        hdf_dims_stream << num_timesteps << " " << num_nodes << " " << this->mNumVariables;
        p_hdf_element->setAttribute(X("Dimensions"), X(hdf_dims_stream.str()));
        p_hype_element->appendChild(p_hdf_element);

        /*
         * e.g. ../cube_2mm_12_elements.h5:/Data
         */
        std::stringstream HDFStream;
        HDFStream << "../" << AbstractHdf5Converter<ELEMENT_DIM, SPACE_DIM>::mFileBaseName << ".h5:/Data";
        DOMText* p_hdf_text = pDomDocument->createTextNode(X(HDFStream.str()));
        p_hdf_element->appendChild(p_hdf_text);
    }
}

#endif

// Explicit instantiation
template class Hdf5ToXdmfConverter<1,1>;
template class Hdf5ToXdmfConverter<1,2>;
template class Hdf5ToXdmfConverter<2,2>;
template class Hdf5ToXdmfConverter<1,3>;
template class Hdf5ToXdmfConverter<2,3>;
template class Hdf5ToXdmfConverter<3,3>;
