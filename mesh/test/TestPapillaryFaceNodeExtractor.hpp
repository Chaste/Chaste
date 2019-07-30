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


#ifndef TESTPAPILLARYFACENODEEXTRACTOR_HPP_
#define TESTPAPILLARYFACENODEEXTRACTOR_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TrianglesMeshReader.hpp"
#include "TetrahedralMesh.hpp"

#include "FakePetscSetup.hpp"

class TestPapillaryExtractor : public CxxTest::TestSuite
{
public:

    void TestExtractor()
    {
//        TrianglesMeshReader<3, 3> mesh_reader("/home/chaste/heart_data/heartT_renum_i");
//        TetrahedralMesh<3,3> mesh;
//        mesh.ConstructFromMeshReader(mesh_reader);
//        TS_ASSERT_EQUALS(mesh.GetNumElements(), 24217344u);
        //TrianglesMeshReader<3, 3> mesh_reader("notforrelease/test/data/mbishop_rabbit_mesh/epi_n_hr.tri");

        //TrianglesMeshReader<3, 3> mesh_reader("notforrelease/test/data/mbishop_rabbit_mesh/rv_n_hr.tri");

        std::string epi_face_file = "/home/chaste/heart_data/pap_face_n.tri";

        std::ifstream nodesfile;
        std::ifstream coordsfile;
        std::ifstream pap_facefile;
        nodesfile.open("/home/chaste/heart_data/heartT_renum_i.tetras");
        coordsfile.open("/home/chaste/heart_data/heartT_renum_i.pts");
        pap_facefile.open("/home/chaste/heart_data/pap_face_n.dat");

        /////////////////////////////////////////////////////////////
        // Defines the numbers of nodes and face nodes
        /////////////////////////////////////////////////////////////
        int num_nodes = 4310704;
        int num_elements = 24217344;
        int num_pap_face = 42539;

        // Defines element list
        int **nodes;
        nodes = new int *[num_elements];
        for (int i=0;i<num_elements;i++)
        {
            nodes[i] = new int[5];
        }
        // Defines coordinate list
        double **coords;
        coords = new double *[num_nodes];
        for (int i=0;i<num_nodes;i++)
        {
            coords[i] = new double[3];
        }
        // Defines list of all papillary face nodes
        int* pap_face;
        pap_face = new int [num_pap_face];


        /////////////////////////////////////////////////////////////////
        // Reads the mesh details from files into matrices and vectors which can be used
        /////////////////////////////////////////////////////////////////
        // Reads in elements
        int node0, node1, node2, node3, tag, dummy;
        nodesfile >> dummy;
        for (int i=0;i<num_elements;i++)
        {
            nodesfile >> node0 >> node1 >> node2 >> node3 >> tag;
            nodes[i][0] = node0;
            nodes[i][1] = node1;
            nodes[i][2] = node2;
            nodes[i][3] = node3;
            nodes[i][4] = tag;
        }
        std::cout << nodes[10][0] << " " << nodes[10][4] << std::endl;

        // Defines scaling factors for coords file
        double x_factor = 0.053/1000;
        double y_factor = x_factor;
        double z_factor = 0.049/1000;

        // Reads in coords file
        double x,y,z;
        coordsfile >> dummy;
        for (int i=0;i<num_nodes;i++)
        {
            coordsfile >> x >> y >> z;
            coords[i][0] = x*x_factor;
            coords[i][1] = y*y_factor;
            coords[i][2] = z*z_factor;

        }
        std::cout << coords[10][0] << " " << coords[10][2] << "\n ";


        // Reads in list of papillary nodes
        int pap_face_value;
        for (int i=0;i<num_pap_face;i++)
        {
            pap_facefile >> pap_face_value;
            pap_face[i] = pap_face_value;
        }
        std::cout << pap_face[10] << "\n ";

        // Defines a list into which radial vectors are stored
        double **gradients;
        gradients = new double *[num_elements];
        for (int i=0;i<num_elements;i++)
        {
            gradients[i] = new double[3];
        }

        // Defines quantities used below
        double x_c,y_c,z_c,x_f,y_f,z_f;
        double r,r_temp;
        int nearest_face = 0;


        // Loops over all elements finding radius vector
        for (int i=0;i<num_elements;i++)
        {
            // Checks to see if we're in a papillary element
            if (nodes[i][4] == 6)
            {
                // Defines coordinates of the centroid of the element
                x_c = 0.25*(coords[nodes[i][0]][0] + coords[nodes[i][1]][0] + coords[nodes[i][2]][0] + coords[nodes[i][3]][0]);
                y_c = 0.25*(coords[nodes[i][0]][1] + coords[nodes[i][1]][1] + coords[nodes[i][2]][1] + coords[nodes[i][3]][1]);
                z_c = 0.25*(coords[nodes[i][0]][2] + coords[nodes[i][1]][2] + coords[nodes[i][2]][2] + coords[nodes[i][3]][2]);

                // Sets the distance to be very big by default
                r = 1000000;

                // Loops over all papillary face nodes
                for (int j=0;j<num_pap_face;j++)
                {
                    // Defines the coordinates of the papillary face node
                    x_f = coords[pap_face[j]][0];
                    y_f = coords[pap_face[j]][1];
                    z_f = coords[pap_face[j]][2];

                    // Calculates the distance between the papillary face node and the centroid
                    r_temp = sqrt( (x_c-x_f)*(x_c-x_f) + (y_c-y_f)*(y_c-y_f) + (z_c-z_f)*(z_c-z_f));

                    // Checks to see if it is the smallest so far - if it is, update the current smallest distance
                    if (r_temp < r)
                    {
                        r = r_temp;
                        nearest_face = j;
                    }

                }
                // Once we have the papillary face node which is closest, use this to re-define its coordinates
                x_f = coords[pap_face[nearest_face]][0];
                y_f = coords[pap_face[nearest_face]][1];
                z_f = coords[pap_face[nearest_face]][2];

                // Defines the radius vector to be r = r1 - r2
                gradients[i][0] = x_c - x_f;
                gradients[i][1] = y_c - y_f;
                gradients[i][2] = z_c - z_f;

            }
        }


        // Writes-out the radius vector file
        std::ofstream vectorfile("/home/chaste/heart_data/radius_vector.dat");
        for (int i=0;i<num_elements;i++)
        {
            vectorfile << gradients[i][0] << " " << gradients[i][1] << " " << gradients[i][2] << "\n";
        }
        vectorfile.close();

        std::cout << "done! \n";

        nodesfile.close();
        coordsfile.close();
        pap_facefile.close();
    }
};

#endif /*TESTPAPILLARYFACENODEEXTRACTOR_HPP_*/
