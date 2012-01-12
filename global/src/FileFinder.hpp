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

#ifndef FILEFINDER_HPP_
#define FILEFINDER_HPP_

#include <string>

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
