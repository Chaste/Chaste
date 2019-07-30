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

#ifndef ABSTRACTDYNAMICALLYLOADABLEENTITY_HPP_
#define ABSTRACTDYNAMICALLYLOADABLEENTITY_HPP_

#include "DynamicCellModelLoader.hpp"
#include "DynamicModelLoaderRegistry.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/split_member.hpp>

/**
 * A mixin class for things that get loaded dynamically to maintain an association between instance objects
 * and the shared library they have been loaded from.
 */
class AbstractDynamicallyLoadableEntity
{
public:
    /**
     * @return a shared pointer to the loader
     */
    const DynamicCellModelLoaderPtr GetLoader() const;

    /**
     * Should only be called by a dynamic cell model loader
     *
     * @param pLoader a shared pointer to the loader
     */
    void SetLoader(DynamicCellModelLoaderPtr pLoader);

    /** Virtual destructor to ensure we're polymorphic */
    virtual ~AbstractDynamicallyLoadableEntity();

private:

    /** The loader for our shared object file */
    DynamicCellModelLoaderPtr mpLoader;

// LCOV_EXCL_START
    /* The save and load methods in this mixin class are both adequately covered in TestDynamicallyLoadedCellModels::TestArchiving()
     * and also in higher-level TestCardiacSimulation tests.  This coverage does not appear in standard gcov parsing and
     * hence spuriously fails.
     */

    friend class boost::serialization::access;
    /**
     * Save the path to the loadable module.
     *
     * @param archive the archive
     * @param version the archive version
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        const std::string so_path = GetLoader()->GetLoadableModulePath();
        archive & so_path;
    }
    /**
     * Load the path to the loadable module, and set our loader from the registry.
     *
     * @param archive the archive
     * @param version the archive version
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        std::string so_path;
        archive & so_path;
        SetLoader(DynamicModelLoaderRegistry::Instance()->GetLoader(so_path));
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
// LCOV_EXCL_STOP
};

#endif /*ABSTRACTDYNAMICALLYLOADABLEENTITY_HPP_*/


