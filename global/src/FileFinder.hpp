/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef FILEFINDER_HPP_
#define FILEFINDER_HPP_

#include <string>

#include "BoostFilesystem.hpp"

/**
 * Structure encapsulating the enumeration of path 'types', i.e. what a path
 * can be relative to.  This allows us to write things like RelativeTo::ChasteTestOutput
 * for readability.
 */
struct RelativeTo
{
    /**
     * What things a path can be relative to.
     */
    enum Value
    {
        CWD, /**< The current working directory */
        ChasteTestOutput, /**< The CHASTE_TEST_OUTPUT folder */
        ChasteSourceRoot, /**< The root of the Chaste source tree */
        Absolute, /**< This is an absolute path */
        AbsoluteOrCwd /**< If it starts with a / then it's absolute, otherwise interpret relative to CWD */
    };
};

/**
 * A helper class for finding files or directories, given paths which can be relative
 * to various locations (e.g. the Chaste source tree root, the current directory, the
 * Chaste test output directory, or an absolute path).
 */
class FileFinder
{
private:

    /** The absolute path to our file. */
    std::string mAbsPath;

    /** Whether to fake one of the fixed paths, e.g. ChasteSourceRoot. */
    static bool msFaking;

    /** Which path to fake. */
    static RelativeTo::Value msFakeWhat;

    /** The fake value of the faked path. */
    static std::string msFakePath;

public:

    /**
     * Default constructor for subclasses to use. They @b must call
     * SetAbsolutePath() in their constructor.
     *
     * This also allows classes to store a FileFinder instance that hasn't
     * been properly set up yet, and assign to it later using operator=.
     */
    FileFinder();

    /**
     * Main constructor.
     *
     * @param rPath  the path to the file/dir to find
     * @param relativeTo  how to interpret this path
     */
    FileFinder(const std::string& rPath, RelativeTo::Value relativeTo);

    /**
     * Find a file (or folder) relative to some file or directory.
     * If the second argument is a directory, we look for the given leaf name within it.
     * If the second argument is a file, then we look for a sibling.
     * An exception is raised if rParentOrSibling does not exist.
     *
     * @param rLeafName  the leaf name of the file/dir to find
     * @param rParentOrSibling  where to look for it
     */
    FileFinder(const std::string& rLeafName, const FileFinder& rParentOrSibling);

    /**
     * Conversion constructor from a Boost Filesystem path object.
     * Note that since fs::path has a conversion constructor from std::string,
     * this allows us to be initialised with a string or character constant, too.
     * The path will be interpreted as relative to the current working directory,
     * unless it is an absolute path.
     *
     * @param rPath  the path to the file/dir to find
     */
    FileFinder(const fs::path& rPath);

    /**
     * Change this FileFinder to point at a new location.
     *
     * @param rPath  the path to the file/dir to find
     * @param relativeTo  how to interpret this path
     */
    void SetPath(const std::string& rPath, RelativeTo::Value relativeTo);

    /**
     * Change this FileFinder to point at a new location, relative to some file or directory.
     *
     * @param rLeafName  the leaf name of the file/dir to find
     * @param rParentOrSibling  where to look for it
     */
    void SetPath(const std::string& rLeafName, const FileFinder& rParentOrSibling);

    /**
     * Test whether we exist.
     */
    bool Exists() const;

    /**
     * Are we pointing at a file?
     */
    bool IsFile() const;

    /**
     * Are we pointing at a directory?
     */
    bool IsDir() const;

    /**
     * Get the absolute path to this file/dir.
     *
     * If this is a directory that exists, the absolute path is guaranteed to end in a '/'.
     * If the directory doesn't exist, it will depend on what was supplied to the constructor.
     */
    std::string GetAbsolutePath() const;

    /**
     * Test whether this file/dir is newer than another file/dir.
     * Compares modification times.
     *
     * @param rOtherEntity  the entity to test against.
     */
    bool IsNewerThan(const FileFinder& rOtherEntity) const;

    /**
     * Get the leaf name of this file or directory.
     *
     * i.e. the individual file or directory name and none of the preceding folders on its path.
     */
    std::string GetLeafName() const;

    /**
     * Get the leaf name of this file or directory, with any file extension removed.
     *
     * i.e. the individual file or directory name and none of the preceding folders on its path.
     */
    std::string GetLeafNameNoExtension() const;

    /**
     * Get a finder for the folder containing this file or directory.
     */
    FileFinder GetParent() const;

    /**
     * Recursively remove this file or folder.
     * Since this is a potentially very dangerous operation, only locations under the Chaste
     * test output folder may be removed.
     */
    void Remove() const;

    /**
     * Test whether a path is absolute. Currently just checks whether the first character is '/'.
     *
     * @param rPath The path to test
     */
    static bool IsAbsolutePath(const std::string& rPath);

    /**
     * Replace any spaces in a path or filename with underscores.
     *
     * @param rPath  A path or file name
     */
    static void ReplaceSpacesWithUnderscores(std::string& rPath);

    /**
     * Replace any underscores in a path or filename with spaces (for making titles etc.).
     *
     * @param rPath  A path or file name
     */
    static void ReplaceUnderscoresWithSpaces(std::string& rPath);

    /**
     * For testing purposes, fake the value of one of the normally fixed paths, e.g. ChasteSourceRoot.
     *
     * @param fakeWhat  which path to fake
     * @param rFakePath  its fake value
     */
    static void FakePath(RelativeTo::Value fakeWhat, const std::string& rFakePath);

    /**
     * Stop faking one of the fixed paths.
     */
    static void StopFaking();
};

#endif /*FILEFINDER_HPP_*/
