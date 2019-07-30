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

#include <sstream>
#include <map>

#include "XdmfMeshWriter.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "Version.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
XdmfMeshWriter<ELEMENT_DIM, SPACE_DIM>::XdmfMeshWriter(const std::string& rDirectory,
                                                       const std::string& rBaseName,
                                                       const bool clearOutputDir)
    : AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory, rBaseName, clearOutputDir),
      mNumberOfTimePoints(1u),
      mTimeStep(1.0)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void XdmfMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMesh(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                                                                 bool keepOriginalElementIndexing)
{
#ifdef _MSC_VER
    EXCEPTION("XDMF is not supported under Windows at present.");
#else
    assert(keepOriginalElementIndexing);
    this->mpDistributedMesh = dynamic_cast<DistributedTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* >(&rMesh);
    bool mesh_is_distributed = (this->mpDistributedMesh != nullptr) && PetscTools::IsParallel();

    if (PetscTools::AmMaster())
    {
        // Write main test Grid collection (to be later replaced by temporal collection)
        // Write references to geometry and topology chunk(s)
        unsigned num_chunks = 1;
        if (mesh_is_distributed)
        {
            num_chunks = PetscTools::GetNumProcs();
        }
        WriteXdmfMasterFile(num_chunks);
    }
    if (!mesh_is_distributed && !PetscTools::AmMaster())
    {
        //If the mesh is not distributed then the master knows everything and will write the geometry/topology as a single chunk
        PetscTools::Barrier("XdmfMeshWriter wait for chunks to be written");
        return;
    }

    // Geometry
    std::stringstream local_geometry_file_name;
    local_geometry_file_name << this->mBaseName << "_geometry_"<< PetscTools::GetMyRank() <<".xml";
    out_stream geometry_file = this->mpOutputFileHandler->OpenOutputFile(local_geometry_file_name.str());
    std::string geom_type = "XYZ";
    if (SPACE_DIM == 2)
    {
        geom_type = "XY";
    }
    (*geometry_file) << "<Geometry GeometryType=\""<< geom_type <<"\">\n";
    unsigned num_nodes = rMesh.GetNumNodes();
    if (this->mpDistributedMesh)
    {
        num_nodes = this->mpDistributedMesh->GetNumLocalNodes() + this->mpDistributedMesh->GetNumHaloNodes();
    }

    (*geometry_file) << "\t<DataItem Format=\"XML\" Dimensions=\""<< num_nodes <<" "<< SPACE_DIM <<"\" DataType=\"Float\">";

    //  Map a global node index into a local index (into mNodes and mHaloNodes as if they were concatenated)
    std::map<unsigned, unsigned> global_to_node_index_map;

    //Node index that we are writing to the chunk (index into mNodes and mHaloNodes as if they were concatenated)
    unsigned index = 0;

    // Owned nodes come first
    for (typename AbstractMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator iter = rMesh.GetNodeIteratorBegin();
         iter != rMesh.GetNodeIteratorEnd();
         ++iter)
    {
        global_to_node_index_map[iter->GetIndex()] = index;
        index++;
        (*geometry_file) << "\n\t\t";
        c_vector<double, SPACE_DIM> current_item = (iter)->rGetLocation();
        for (unsigned j=0; j<SPACE_DIM; j++)
        {
            (*geometry_file) << current_item[j] << "\t";
        }
    }

    // Halo nodes
    if (this->mpDistributedMesh)
    {
        for (typename DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::HaloNodeIterator halo_iter=this->mpDistributedMesh->GetHaloNodeIteratorBegin();
                halo_iter != this->mpDistributedMesh->GetHaloNodeIteratorEnd();
                ++halo_iter)
        {
            global_to_node_index_map[(*halo_iter)->GetIndex()] = index;
            index++;
            (*geometry_file) << "\n\t\t";
            c_vector<double, SPACE_DIM> current_item = (*halo_iter)->rGetLocation();
            for (unsigned j=0; j<SPACE_DIM; j++)
            {
                (*geometry_file) << current_item[j] << "\t";
            }
        }
    }
    (*geometry_file) << "\n";

    (*geometry_file) << "\t</DataItem>\n";
    (*geometry_file) << "</Geometry>\n";
    (*geometry_file) << "<!-- " + ChasteBuildInfo::GetProvenanceString() + "-->\n";
    geometry_file->close();

    // Topology
    std::stringstream local_topology_file_name;
    local_topology_file_name << this->mBaseName << "_topology_"<< PetscTools::GetMyRank() <<".xml";
    out_stream topology_file = this->mpOutputFileHandler->OpenOutputFile(local_topology_file_name.str());
    std::string top_type = "Tetrahedron";
    if (SPACE_DIM == 2)
    {
        top_type = "Triangle";
    }
    unsigned num_elems = rMesh.GetNumElements();
    if (this->mpDistributedMesh)
    {
        num_elems = this->mpDistributedMesh->GetNumLocalElements();
    }
    (*topology_file) << "<Topology TopologyType=\""<< top_type <<"\" NumberOfElements=\""<< num_elems <<"\">\n";
    (*topology_file) << "\t<DataItem Format=\"XML\" Dimensions=\""<< num_elems <<" "<< ELEMENT_DIM+1 <<"\">";
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator elem_iter = rMesh.GetElementIteratorBegin();
         elem_iter != rMesh.GetElementIteratorEnd();
         ++elem_iter)
    {
        (*topology_file) << "\n\t\t";
        for (unsigned j=0; j<ELEMENT_DIM+1; j++)
        {
            unsigned local_index = global_to_node_index_map[ elem_iter->GetNodeGlobalIndex(j) ];
            (*topology_file) << local_index <<"\t";
        }
    }
    (*topology_file) << "\n";

    (*topology_file) << "\t</DataItem>\n";
    (*topology_file) << "</Topology>\n";
    (*topology_file) << "<!-- " + ChasteBuildInfo::GetProvenanceString() + "-->\n";
    topology_file->close();
    PetscTools::Barrier("XdmfMeshWriter wait for chunks to be written");
#endif
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void XdmfMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFiles()
{
#ifdef _MSC_VER
    EXCEPTION("XDMF is not supported under Windows at present.");
#else
    // This method is only called when there is no mesh.  We are writing from a reader.
    if (PetscTools::AmMaster())
    {
        WriteXdmfMasterFile();

        // Geometry
        out_stream geometry_file = this->mpOutputFileHandler->OpenOutputFile(this->mBaseName + "_geometry_0.xml");
        std::string geom_type = "XYZ";
        if (SPACE_DIM == 2)
        {
            geom_type = "XY";
        }
        (*geometry_file) << "<Geometry GeometryType=\""<< geom_type <<"\">\n";
        (*geometry_file) << "\t<DataItem Format=\"XML\" Dimensions=\""<< this->GetNumNodes() <<" "<< SPACE_DIM <<"\" DataType=\"Float\">";
        for (unsigned item_num=0; item_num<this->GetNumNodes(); item_num++)
        {
            (*geometry_file) << "\n\t\t";
            std::vector<double> current_item = this->GetNextNode();
            for (unsigned j=0; j<SPACE_DIM; j++)
            {
                (*geometry_file) << current_item[j]<<"\t";
            }
        }
        (*geometry_file) << "\n";

        (*geometry_file) << "\t</DataItem>\n";
        (*geometry_file) << "</Geometry>\n";
        (*geometry_file) << "<!-- " + ChasteBuildInfo::GetProvenanceString() + "-->\n";
        geometry_file->close();

        // Topology
        out_stream topology_file = this->mpOutputFileHandler->OpenOutputFile(this->mBaseName + "_topology_0.xml");
        std::string top_type = "Tetrahedron";
        if (SPACE_DIM == 2)
        {
            top_type = "Triangle";
        }
        (*topology_file) << "<Topology TopologyType=\""<< top_type <<"\" NumberOfElements=\""<< this->GetNumElements() <<"\">\n";
        (*topology_file) << "\t<DataItem Format=\"XML\" Dimensions=\""<< this->GetNumElements() <<" "<< ELEMENT_DIM+1 <<"\">";
        for (unsigned item_num=0; item_num<this->GetNumElements(); item_num++)
        {
            (*topology_file) << "\n\t\t";
            std::vector<unsigned> current_item = this->GetNextElement().NodeIndices;
            for (unsigned j=0; j<ELEMENT_DIM+1; j++)
            {
                (*topology_file) << current_item[j]<<"\t";
            }
        }
        (*topology_file) << "\n";

        (*topology_file) << "\t</DataItem>\n";
        (*topology_file) << "</Topology>\n";
        (*topology_file) << "<!-- " + ChasteBuildInfo::GetProvenanceString() + "-->\n";
        topology_file->close();
    }
#endif
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void XdmfMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteXdmfMasterFile(unsigned numberOfChunks)
{
#ifndef _MSC_VER
    assert(PetscTools::AmMaster());
    // Define namespace symbols
    XERCES_CPP_NAMESPACE_USE

    // Initialize Xerces.
    XMLPlatformUtils::Initialize();

    DOMImplementation* p_DOM_implementation = DOMImplementationRegistry::getDOMImplementation(X("core"));

    DOMDocumentType* p_DOM_document_type = p_DOM_implementation->createDocumentType(X("Xdmf"),nullptr,X("Xdmf.dtd"));
    DOMDocument* p_DOM_document = p_DOM_implementation->createDocument(nullptr, X("Xdmf"), p_DOM_document_type);
    DOMElement* p_root_element = p_DOM_document->getDocumentElement();
    p_root_element->setAttribute(X("Version"), X("2.0"));
    p_root_element->setAttribute(X("xmlns:xi"), X("http://www.w3.org/2001/XInclude"));

    DOMElement* p_domain_element =  p_DOM_document->createElement(X("Domain"));
    p_root_element->appendChild(p_domain_element);

    // Temporal collection
    DOMElement* p_grid_temp_collection_element =  p_DOM_document->createElement(X("Grid"));
    p_grid_temp_collection_element->setAttribute(X("CollectionType"), X("Temporal"));
    p_grid_temp_collection_element->setAttribute(X("GridType"), X("Collection"));
    p_domain_element->appendChild(p_grid_temp_collection_element);

    // Time values
    DOMElement* p_time_element =  p_DOM_document->createElement(X("Time"));
    p_time_element->setAttribute(X("TimeType"), X("HyperSlab"));
    p_grid_temp_collection_element->appendChild(p_time_element);

    DOMElement* p_time_dataitem_element =  p_DOM_document->createElement(X("DataItem"));
    p_time_dataitem_element->setAttribute(X("Format"),X("XML"));
    p_time_dataitem_element->setAttribute(X("NumberType"),X("Float"));
    p_time_dataitem_element->setAttribute(X("Dimensions"),X("3"));
    p_time_element->appendChild(p_time_dataitem_element);

    std::stringstream time_stream;
    time_stream << "0.0 " << mTimeStep << " " << mNumberOfTimePoints;
    DOMText* p_time_text = p_DOM_document->createTextNode(X(time_stream.str()));
    p_time_dataitem_element->appendChild(p_time_text);

    for (unsigned t=0; t<mNumberOfTimePoints; ++t)
    {
        DOMElement* p_grid_collection_element =  p_DOM_document->createElement(X("Grid"));
        p_grid_collection_element->setAttribute(X("CollectionType"), X("Spatial"));
        p_grid_collection_element->setAttribute(X("GridType"), X("Collection"));
        //p_grid_collection_element->setAttribute(X("Name"), X("spatial_collection"));
        p_grid_temp_collection_element->appendChild(p_grid_collection_element);

        if (t==0)
        {
            for (unsigned chunk=0; chunk<numberOfChunks; chunk++)
            {
                std::stringstream chunk_stream;
                chunk_stream << chunk;

                DOMElement* p_grid_element =  p_DOM_document->createElement(X("Grid"));
                p_grid_element->setAttribute(X("GridType"), X("Uniform"));
                p_grid_element->setAttribute(X("Name"), X("Chunk_" + chunk_stream.str()));
                p_grid_collection_element->appendChild(p_grid_element);

                //DOMElement* p_geom_element =  p_DOM_document->createElement(X("Geometry"));
                //p_geom_element->setAttribute(X("Reference"),X("/Xdmf/Domain/Geometry[1]"));
                DOMElement* p_geom_element =  p_DOM_document->createElement(X("xi:include"));
                p_geom_element->setAttribute(X("href"), X(this->mBaseName+"_geometry_"+chunk_stream.str()+".xml"));
                p_grid_element->appendChild(p_geom_element);
                //DOMElement* p_topo_element =  p_DOM_document->createElement(X("Topology"));
                //p_topo_element->setAttribute(X("Reference"),X("/Xdmf/Domain/Topology[1]"));
                DOMElement* p_topo_element =  p_DOM_document->createElement(X("xi:include"));
                p_topo_element->setAttribute(X("href"), X(this->mBaseName+"_topology_"+chunk_stream.str()+".xml"));
                p_grid_element->appendChild(p_topo_element);

                /*
                 * p_grid_element may now need an Attribute (node data). Call Annotate,
                 * which here does nothing, but in pde can be overloaded to print variables
                 */
                AddDataOnNodes(p_grid_element, p_DOM_document, t);
            }
        }
        else // t>0
        {
            for (unsigned chunk=0; chunk<numberOfChunks; chunk++)
            {
                std::stringstream chunk_stream;
                chunk_stream << chunk;

                DOMElement* p_grid_element =  p_DOM_document->createElement(X("Grid"));
                p_grid_element->setAttribute(X("GridType"), X("Subset"));
                p_grid_element->setAttribute(X("Section"), X("All"));
                p_grid_collection_element->appendChild(p_grid_element);

                /*
                 * p_grid_element may now need an Attribute (node data). Call Annotate,
                 * which here does nothing, but in pde can be overloaded to print variables
                 */
                AddDataOnNodes(p_grid_element, p_DOM_document, t);
                DOMElement* p_grid_ref_element =  p_DOM_document->createElement(X("Grid"));
                p_grid_ref_element->setAttribute(X("GridType"), X("Uniform"));
                p_grid_ref_element->setAttribute(X("Reference"), X("XML"));
                //p_grid_ref_element->setAttribute(X("Name"), X("Chunk_" + chunk_stream.str()));

                DOMText* p_ref_text = p_DOM_document->createTextNode(X("/Xdmf/Domain/Grid/Grid/Grid[@Name=\"Chunk_"+chunk_stream.str()+"\"]"));
                p_grid_ref_element->appendChild(p_ref_text);
                p_grid_element->appendChild(p_grid_ref_element);
            }
        }
    }
    // Create a Comment node, and then append this to the root element.
    DOMComment* p_provenance_comment = p_DOM_document->createComment(X(" "+ChasteBuildInfo::GetProvenanceString()));
    p_DOM_document->appendChild(p_provenance_comment);

    XMLFormatTarget* p_target = new LocalFileFormatTarget(X(this->mpOutputFileHandler->GetOutputDirectoryFullPath() + this->mBaseName+".xdmf"));

#if _XERCES_VERSION >= 30000
    DOMLSSerializer* p_serializer = ((DOMImplementationLS*)p_DOM_implementation)->createLSSerializer();
    p_serializer->getDomConfig()->setParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true);
    DOMLSOutput* p_output = ((DOMImplementationLS*)p_DOM_implementation)->createLSOutput(); // Calls a new somewhere!
    p_output->setByteStream(p_target);
    p_serializer->write(p_DOM_document, p_output);
#else
    DOMWriter* p_serializer = ((DOMImplementationLS*)p_DOM_implementation)->createDOMWriter();
    p_serializer->setFeature(XMLUni::fgDOMWRTFormatPrettyPrint, true);
    p_serializer->writeNode(p_target, *p_DOM_document);
#endif

    // Cleanup
    p_serializer->release();
    p_DOM_document->release();
#if _XERCES_VERSION >= 30000
    delete p_output;
#endif
    delete p_target;
    XMLPlatformUtils::Terminate();
#endif // _MSC_VER
}

// Explicit instantiation
template class XdmfMeshWriter<1,1>;
template class XdmfMeshWriter<1,2>;
template class XdmfMeshWriter<1,3>;
template class XdmfMeshWriter<2,2>; // Actually used
template class XdmfMeshWriter<2,3>;
template class XdmfMeshWriter<3,3>; // Actually used
