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

#ifndef TESTRANDOMNUMBERGENERATOR_HPP_
#define TESTRANDOMNUMBERGENERATOR_HPP_
#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "OutputFileHandler.hpp"
#include "RandomNumberGenerator.hpp"
// We use a c_matrix for a bit of storage.  It's rather naughty using a linalg header
// here, but since it doesn't have a corresponding .cpp file we get away with it!
#include "UblasMatrixInclude.hpp"

class TestRandomNumberGenerator : public CxxTest::TestSuite
{
public:

    double ran1;

    void TestRandomNumbers()
    {
        srandom(0);
        ran1 = (double)random()/RAND_MAX;

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

        double ran2 = p_gen->ranf();
        TS_ASSERT_DELTA(ran1, ran2, 1e-7);

        RandomNumberGenerator::Destroy();
    }

    void TestNewMethodSeed()
    {
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        double ran2 = p_gen->ranf();
        TS_ASSERT_DELTA(ran1, ran2, 1e-7);

        RandomNumberGenerator::Destroy();
    }


    void TestOtherRandomStuffDestroysRandomSequence()
    {
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

        //First reseed and get the first random number (as above)
        p_gen->Reseed(0);
        double ran2 = p_gen->ranf();
        TS_ASSERT_DELTA(ran1, ran2, 1e-7);

        //Now reseed, do something else random and then get "the first" random number
        p_gen->Reseed(0);
        std::vector<unsigned> some_vector(10);
        std::random_shuffle(some_vector.begin(), some_vector.end());
        double ran3 = p_gen->ranf();
        TS_ASSERT_DIFFERS(ran1, ran3);

        //Again - with rand()
        p_gen->Reseed(0);
        rand();
        double ran4 = p_gen->ranf();
        TS_ASSERT_DIFFERS(ran1, ran4);

        //Again - with random()
        p_gen->Reseed(0);
        random();
        double ran5 = p_gen->ranf();
        TS_ASSERT_DIFFERS(ran1, ran5);

        //Again - with nothing
        p_gen->Reseed(0);
        double ran6 = p_gen->ranf();
        TS_ASSERT_DELTA(ran1, ran6, 1e-7);

        RandomNumberGenerator::Destroy();

    }

    void TestDifferentRandomSeed()
    {
        srandom(36);
        ran1 = (double)random()/RAND_MAX;

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        p_gen->Reseed(36);


        double ran2 = p_gen->ranf();
        TS_ASSERT_DELTA(ran1, ran2, 1e-7);

        RandomNumberGenerator::Destroy();
    }

    void TestArchiveRandomNumberGenerator()
    {
        OutputFileHandler handler("archive",false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "random_number.arch";

        std::vector<double> generated_numbers;

        // Create and archive random number generator
        {
            // Save random number generator
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            p_gen->Reseed(7); // This gives us full coverage

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            for (unsigned i=0; i<5; i++)
            {
                p_gen->ranf();
            }

            // A few extra calls before archiving
            p_gen->ranf();
            p_gen->randMod(3);
            p_gen->StandardNormalRandomDeviate();
            p_gen->NormalRandomDeviate(0.5, 0.1);

            SerializableSingleton<RandomNumberGenerator>* const p_wrapper = p_gen->GetSerializationWrapper();
            output_arch << p_wrapper;

            // Make sure saving it twice gets the same instance
            {
                SerializableSingleton<RandomNumberGenerator>* const p_wrapper = p_gen->GetSerializationWrapper();
                output_arch << p_wrapper;
            }

            // Generator saved here - record the next 10 numbers
            for (unsigned i=0; i<10; i++)
            {
                double random = p_gen->ranf();
                generated_numbers.push_back(random);
            }

            /*
             * Rcord some numbers from the normal distribution too.
             * We generate quite a few for coverage of the three cases
             * in RandomNumberGenerator::StandardNormalRandomDeviate().
             */
            for (unsigned i=0; i<10; i++)
            {
                double random = p_gen->NormalRandomDeviate(0.5, 0.1);
                generated_numbers.push_back(random);
            }

            RandomNumberGenerator::Destroy();
        }

        // Restore
        {
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            p_gen->Reseed(25); // any old seed
            for (unsigned i=0; i<7; i++) // generate some numbers
            {
                p_gen->ranf();
            }

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            SerializableSingleton<RandomNumberGenerator>* p_orig_wrapper = p_gen->GetSerializationWrapper();
            {
                SerializableSingleton<RandomNumberGenerator>* p_wrapper;
                input_arch >> p_wrapper;
                TS_ASSERT_DIFFERS(p_wrapper, p_gen->GetSerializationWrapper());
                TS_ASSERT_EQUALS(p_orig_wrapper, p_gen->GetSerializationWrapper());
                TS_ASSERT_EQUALS(p_gen, RandomNumberGenerator::Instance());
            }

            {
                SerializableSingleton<RandomNumberGenerator>* p_wrapper;
                input_arch >> p_wrapper;
                TS_ASSERT_DIFFERS(p_wrapper, p_gen->GetSerializationWrapper());
                TS_ASSERT_EQUALS(p_orig_wrapper, p_gen->GetSerializationWrapper());
                TS_ASSERT_EQUALS(p_gen, RandomNumberGenerator::Instance());
            }

            /*
             * Random Number generator restored.
             * Check it generates the same numbers as the one we saved.
             */
            for (unsigned i=0; i<generated_numbers.size(); i++)
            {
                double random;
                if (i<10)
                {
                    random = p_gen->ranf();
                }
                else
                {
                    random = p_gen->NormalRandomDeviate(0.5, 0.1);
                }
                TS_ASSERT_DELTA(random,generated_numbers[i],1e-7);
            }

            RandomNumberGenerator::Destroy();
        }
    }

    void TestShuffle() throw(Exception)
    {
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        p_gen->Reseed(0);

        std::vector<unsigned> shuffled_results;
        p_gen->Shuffle(5, shuffled_results);

        for (unsigned i=0; i<5; i++)
        {
            bool found = false;
            for (unsigned j=0; j<shuffled_results.size(); j++)
            {
                if (shuffled_results[j]==i)
                {
                    found = true;
                    break;
                }
            }
            TS_ASSERT_EQUALS(found, true);
        }

        unsigned num_trials = 100000;
        c_matrix<unsigned,5,5> results = zero_matrix<unsigned>(5,5);

        for (unsigned trial=0; trial<num_trials; trial++)
        {
            p_gen->Shuffle(5, shuffled_results);
            for (unsigned i=0; i<5; i++)
            {
                for (unsigned j=0; j<5; j++)
                {
                    if (shuffled_results[j] == i)
                    {
                        results(i,j)++;
                    }
                }
            }
        }
        for (unsigned i=0; i<5; i++)
        {
            for (unsigned j=0; j<5; j++)
            {
                // Probability of i going to position j
                double prob = (double)results(i,j)/num_trials;

                /*
                 * This test could fail with very low probability (just rerun).
                 * We accept 0.19 to 0.21 (note, usually in 0.199 to 0.201 with
                 * a million trials (we use 10^5 trials))
                 */
                TS_ASSERT_DELTA(prob, 0.2, 1e-2);
            }
        }
    }
};

#endif /*TESTRANDOMNUMBERGENERATOR_HPP_*/
