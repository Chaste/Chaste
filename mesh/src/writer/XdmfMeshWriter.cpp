/*

Copyright (c) 2005-2013, University of Oxford.
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
: AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory, rBaseName, clearOutputDir)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void XdmfMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMesh(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
        bool keepOriginalElementIndexing)
{
    assert(keepOriginalElementIndexing==true);
    this->mpDistributedMesh = dynamic_cast<DistributedTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* >(&rMesh);
    bool mesh_is_distributed = (this->mpDistributedMesh != NULL) && PetscTools::IsParallel();

    if (PetscTools::AmMaster())
    {
        std::string master_file_name = this->mBaseName + ".xdmf";
        out_stream master_file = this->mpOutputFileHandler->OpenOutputFile(master_file_name);

        // Write header information
        (*master_file) << "<?xml version=\"1.0\" ?>\n";
        (*master_file) << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        (*master_file) << "<Xdmf Version=\"2.0\" xmlns:xi=\"http://www.w3.org/2001/XInclude\">\n";
        (*master_file) << "\t<Domain>\n";


        // Write main test Grid collection (to be later replaced by temporal collection)
        // Write references to geometry and topology chunk(s)
        unsigned num_chunks = 1;
        if (mesh_is_distributed)
        {
            num_chunks = PetscTools::GetNumProcs();
        }
        (*master_file) << "\t\t<Grid Name=\"Grid\" GridType=\"Collection\" CollectionType=\"Spatial\">\n";
        for (unsigned chunk=0; chunk<num_chunks; chunk++)
        {
            std::stringstream geometry_file_name;
            geometry_file_name << this->mBaseName << "_geometry_"<< chunk <<".xml";
            std::stringstream topology_file_name;
            topology_file_name << this->mBaseName << "_topology_"<< chunk <<".xml";
            (*master_file) << "\t\t\t<Grid Name=\"Chunk_"<< chunk <<"\" Gridtype=\"Uniform\">\n";
            (*master_file) << "\t\t\t\t<xi:include href=\""<< geometry_file_name.str() <<"\"/>\n";
            (*master_file) << "\t\t\t\t<xi:include href=\""<< topology_file_name.str() <<"\"/>\n";
            (*master_file) << "\t\t\t</Grid>\n";
        }
        (*master_file) << "\t\t</Grid>\n";

        // Write footer
        (*master_file) << "\t</Domain>\n";
        (*master_file) << "</Xdmf>\n";

        (*master_file) << "<!-- " + ChasteBuildInfo::GetProvenanceString() + "-->\n";
        master_file->close();
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
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void XdmfMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFiles()
{
    if (PetscTools::AmMaster())
    {
        std::string master_file_name = this->mBaseName + ".xdmf";
        out_stream master_file = this->mpOutputFileHandler->OpenOutputFile(master_file_name);

        // Write header information
        (*master_file) << "<?xml version=\"1.0\" ?>\n";
        (*master_file) << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        (*master_file) << "<Xdmf Version=\"2.0\" xmlns:xi=\"http://www.w3.org/2001/XInclude\">\n";
        (*master_file) << "\t<Domain>\n";


        // Write main test Grid collection (to be later replaced by temporal collection)
        // Write references to geometry and topology chunk(s)
        std::string geometry_file_name = this->mBaseName + "_geometry_0.xml";
        std::string topology_file_name = this->mBaseName + "_topology_0.xml";
        (*master_file) << "\t\t<Grid Name=\"Chunk_0\" Gridtype=\"Uniform\">\n";
        (*master_file) << "\t\t\t<xi:include href=\""<< geometry_file_name <<"\"/>\n";
        (*master_file) << "\t\t\t<xi:include href=\""<< topology_file_name <<"\"/>\n";
        (*master_file) << "\t\t</Grid>\n";

        // Write footer
        (*master_file) << "\t</Domain>\n";
        (*master_file) << "</Xdmf>\n";
        (*master_file) << "<!-- " + ChasteBuildInfo::GetProvenanceString() + "-->\n";

        master_file->close();

        // Geometry
        out_stream geometry_file = this->mpOutputFileHandler->OpenOutputFile(geometry_file_name);
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
        out_stream topology_file = this->mpOutputFileHandler->OpenOutputFile(topology_file_name);
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
}
/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class XdmfMeshWriter<1,1>;
template class XdmfMeshWriter<1,2>;
template class XdmfMeshWriter<1,3>;
template class XdmfMeshWriter<2,2>; // Actually used
template class XdmfMeshWriter<2,3>;
template class XdmfMeshWriter<3,3>; // Actually used



