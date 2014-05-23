/*

Copyright (c) 2005-2014, University of Oxford.
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

/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 */
#ifndef TESTBIDOMAINWITHCONDUCTIVITYMODIFIERTUTORIAL_HPP_
#define TESTBIDOMAINWITHCONDUCTIVITYMODIFIERTUTORIAL_HPP_
/*
 * = We now show how to run a bidomain simulation using conductivity modifiers. =
 *
 * In this tutorial we run another bidomain simulation,
 * showing (i) an example using one of the source cell factories,
 * (ii) an example of using a conductivity modifier.
 *
 * The first thing to do is to include the headers as before.
 */
#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "LuoRudy1991.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "AbstractConductivityModifier.hpp"
#include "ZeroStimulusCellFactory.hpp"


/*
 * Here we define our conductivity modifier. It inherits from the abstract class.
 * Make sure to use the right element/space dims.
 */
class SimpleConductivityModifier : public AbstractConductivityModifier<2,2>
{

    /*
     * mTensor is "working memory" to hold the returned modified tensor.
     * This is needed because we return by reference to the problem class, so we need
     * to make sure the memory isn't overwritten before it's done with.
     *
     * Most of the time, we will modify the original conductivity tensor,
     * but sometimes we want to return some "constant" matrix.
     */
private:

    c_matrix<double,2,2> mTensor;
    c_matrix<double,2,2> mSpecialMatrix;


public:
     /*
      * Constructor. In Chaste, all conductivity tensors are diagonal, so we initialise our "constant matrix"
      * to zero so that we only need to set the diagonal entries in the constructor.
      */
    SimpleConductivityModifier()
        : AbstractConductivityModifier<2,2>(),
          mSpecialMatrix( zero_matrix<double>(2,2) )
          {
              mSpecialMatrix(0,0) = 3.14;
              mSpecialMatrix(1,1) = 0.707;
          }

   /*
    * `rCalculateModifiedConductivityTensor` returns a reference to the "processed" conductivity tensor.
    */
    c_matrix<double,2,2>& rCalculateModifiedConductivityTensor(unsigned elementIndex, const c_matrix<double,2,2>& rOriginalConductivity, unsigned domainIndex)
    {

        if ( elementIndex == 0 )
        {
            // For element 0 let's return the "special matrix", regardless of intra/extra-cellular.
            return mSpecialMatrix;
        }


        // Otherwise, we change the behaviour depending on the `domainIndex` (intra/extra-cellular).
        double domain_scaling;
        if ( domainIndex == 0 )
        {
            domain_scaling = 1.0; // Intracellular, domainIndex==0
        }
        else
        {
            domain_scaling = 1.5; // Extracellular, domainIndex==1
        }

        // Modify the current conductivity according to some expression by running along the diagonal,
        // save to the "working memory", and return.
        for ( unsigned i=0; i<2; i++ )
        {
            mTensor(i,i) = domain_scaling*elementIndex*rOriginalConductivity(i,i);
        }
        return mTensor;
    }
};


/* Now we define the test class, which must inherit from {{{CxxTest::TestSuite}}}
 * as usual, and the (public) test method
 */
class TestBidomainWithConductivityModifierTutorial : public CxxTest::TestSuite
{
public:
    void TestConductivityModifier() throw(Exception)
    {
        /*
         * Generate the mesh.
         */
        DistributedTetrahedralMesh<2,2> mesh;
        mesh.ConstructRegularSlabMesh(0.5, 1.0, 0.5); // Mesh has 4 elements

        /*
         * Usually you'll provide your own cell factory, here we use a trivial one for simplicity.
         * Then we set up the problem with the factory as usual.
         */
        ZeroStimulusCellFactory<CellLuoRudy1991FromCellML,2> cell_factory;
        BidomainProblem<2> bidomain_problem( &cell_factory );
        bidomain_problem.SetMesh( &mesh );

        /*
         * The problem generates the `BidomainTissue`, but only after `Initialise()`, so do that now and get
         * the tissue so we can apply the modifier later.
         */
        bidomain_problem.Initialise();
        BidomainTissue<2>* p_bidomain_tissue = bidomain_problem.GetBidomainTissue();

        /*
         * Get the original conductivity tensor values. We haven't set them using
         * `HeartConfig->SetIntra/ExtracellularConductivities` so they'll just be the defaults.
         *
         * The first argument below is the element ID (we just check the first element we own here). The second accesses
         * the diagonal elements. We just do (0,0), as (1,1) should be the same (no fibre orientation).
         * Off-diagonal elements will be 0.
         *
         * As we don't have many elements, when we run on more than two processors some processors
         * will not own any elements. We only try to access the conductivity tensors if the process
         * owns at least one element.
         *
         * We then check that we have the correct (default) conductivity values.
         */
        double orig_intra_conductivity_0 = 0.0;
        double orig_extra_conductivity_0 = 0.0;
        if (mesh.GetElementIteratorBegin() != mesh.GetElementIteratorEnd())
        {
            unsigned first_element = mesh.GetElementIteratorBegin()->GetIndex();
            orig_intra_conductivity_0 = p_bidomain_tissue->rGetIntracellularConductivityTensor(first_element)(0,0);
            orig_extra_conductivity_0 = p_bidomain_tissue->rGetExtracellularConductivityTensor(first_element)(0,0);

            TS_ASSERT_DELTA(orig_intra_conductivity_0, 1.75, 1e-9); // hard-coded using default
            TS_ASSERT_DELTA(orig_extra_conductivity_0, 7.0, 1e-9); // hard-coded using default
        }

        /*
         * Now we make the modifier object and set it on the tissue using `SetConductivityModifier`.
         */
        SimpleConductivityModifier modifier;
        p_bidomain_tissue->SetConductivityModifier( &modifier );

        /*
         * Check that the conductivities have changed! We iterate over all elements owned by this process
         * and check their conductivity.
         */
        for (AbstractTetrahedralMesh<2,2>::ElementIterator elt_iter=mesh.GetElementIteratorBegin();
             elt_iter!=mesh.GetElementIteratorEnd();
             ++elt_iter)
        {
            unsigned index = elt_iter->GetIndex();
            if (index == 0u)
            {
                TS_ASSERT_DELTA(p_bidomain_tissue->rGetIntracellularConductivityTensor(0)(0,0), 3.14, 1e-9);
                TS_ASSERT_DELTA(p_bidomain_tissue->rGetExtracellularConductivityTensor(0)(0,0), 3.14, 1e-9);
                TS_ASSERT_DELTA(p_bidomain_tissue->rGetIntracellularConductivityTensor(0)(1,1), 0.707, 1e-9);
                TS_ASSERT_DELTA(p_bidomain_tissue->rGetExtracellularConductivityTensor(0)(1,1), 0.707, 1e-9);
            }
            else
            {
                TS_ASSERT_DELTA(p_bidomain_tissue->rGetIntracellularConductivityTensor(index)(0,0), 1.0*index*orig_intra_conductivity_0, 1e-9);
                TS_ASSERT_DELTA(p_bidomain_tissue->rGetExtracellularConductivityTensor(index)(0,0), 1.5*index*orig_extra_conductivity_0, 1e-9);
            }
        }
    }
};

#endif /*TESTBIDOMAINWITHCONDUCTIVITYMODIFIERTUTORIAL_HPP_*/
