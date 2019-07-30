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


#ifndef SPACECONVERGENCETESTER_HPP_
#define SPACECONVERGENCETESTER_HPP_


#include "AbstractConvergenceTester.hpp"

/**
 * Run the same simulation on cuboid meshes at progressively finer scales
 * until some convergence criterion is met.
 */
template<class CELL, class CARDIAC_PROBLEM, unsigned DIM, unsigned PROBLEM_DIM>
class SpaceConvergenceTester : public AbstractConvergenceTester<CELL, CARDIAC_PROBLEM, DIM, PROBLEM_DIM>
{
public:
    /**
     * The intial mesh is mesh 0 which has a space-step of 0.05cm on a 0.2cm mesh
     */
    void SetInitialConvergenceParameters()
    {
        this->MeshNum=0;
    }
    /**
     * Each new run has an increased the mesh number. This halves the space step and increased
     * the complexity of the simulation by 2^DIM
     */
    void UpdateConvergenceParameters()
    {
        this->MeshNum++;

    }
    /**
     * @return true to give up convergence when the number of unknowns becomes too high (either
     * from memory or CPU perspective).
     */
    bool GiveUpConvergence()
    {
        switch(DIM)
        {
            case 1:
            {
                return this->MeshNum>10;
                break;
            }
            case 2:
            {
                return this->MeshNum>6;
                break;
            }
            case 3:
            {
                return this->MeshNum>4;
                break;
            }
            default:
                NEVER_REACHED;
                return true;//To keep Intel compiler happy
        }
        return true;//To keep Intel compiler happy
    }
    /**
     * @return the typical step size as the abcissa.
     */
    double Abscissa()
    {
        unsigned mesh_size = SmallPow(2u, this->MeshNum+2); // number of elements in each dimension
        return this->mMeshWidth/(double) mesh_size;
    }

    /**
     * @return the number of the current mesh (cast to int).
     */
    int GetMeshNum()
    {
        return (int) this->MeshNum; //unsigned -> int is just cosmetic here.  (The test looks prettier).
    }
    /**
     * @return the space step in the Cartesian directions
     */
    double GetSpaceStep()
    {
        unsigned mesh_size = SmallPow(2u, this->MeshNum+2);// number of elements in each dimension
        double scaling = this->mMeshWidth/(double) mesh_size;
        return scaling;
    }
};

#endif /*SPACECONVERGENCETESTER_HPP_*/
