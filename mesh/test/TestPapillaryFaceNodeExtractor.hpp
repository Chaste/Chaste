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


#ifndef TESTPAPILLARYFACENODEEXTRACTOR_HPP_
#define TESTPAPILLARYFACENODEEXTRACTOR_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TrianglesMeshReader.hpp"
#include "TetrahedralMesh.hpp"
using namespace std;
class TestPapillaryExtractor : public CxxTest::TestSuite
{
public:

    void TestExtractor() throw(Exception)
    {
//        TrianglesMeshReader<3, 3> mesh_reader("/home/chaste/heart_data/heartT_renum_i");
//        TetrahedralMesh<3,3> mesh;
//        mesh.ConstructFromMeshReader(mesh_reader);
//        TS_ASSERT_EQUALS(mesh.GetNumElements(), 24217344u);
        //TrianglesMeshReader<3, 3> mesh_reader("notforrelease/test/data/mbishop_rabbit_mesh/epi_n_hr.tri");

        //TrianglesMeshReader<3, 3> mesh_reader("notforrelease/test/data/mbishop_rabbit_mesh/rv_n_hr.tri");

        std::string epi_face_file = "/home/chaste/heart_data/pap_face_n.tri";

        ifstream nodesfile;
        ifstream coordsfile;
        ifstream pap_facefile;
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
        cout << nodes[10][0] << " " << nodes[10][4] << std::endl;

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
        cout << coords[10][0] << " " << coords[10][2] << "\n ";


        // Reads in list of papillary nodes
        int pap_face_value;
        for (int i=0;i<num_pap_face;i++)
        {
            pap_facefile >> pap_face_value;
            pap_face[i] = pap_face_value;
        }
        cout << pap_face[10] << "\n ";

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
        ofstream vectorfile("/home/chaste/heart_data/radius_vector.dat");
        for (int i=0;i<num_elements;i++)
        {
            vectorfile << gradients[i][0] << " " << gradients[i][1] << " " << gradients[i][2] << "\n";
        }
        vectorfile.close();

        cout << "done! \n";

        nodesfile.close();
        coordsfile.close();
        pap_facefile.close();
    }
};

#endif /*TESTPAPILLARYFACENODEEXTRACTOR_HPP_*/
