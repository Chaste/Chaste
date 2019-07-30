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
        ChasteBuildRoot, /**< The root of the Chaste build tree */
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
public:
    /**
     * Default constructor for subclasses to use. They should call
     * SetPath() in their constructor.
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
     * Needed because we have virtual methods.
     */
    virtual ~FileFinder();

    /**
     * Change this FileFinder to point at a new location.
     *
     * @param rPath  the path to the file/dir to find
     * @param relativeTo  how to interpret this path
     */
    virtual void SetPath(const std::string& rPath, RelativeTo::Value relativeTo);

    /**
     * Change this FileFinder to point at a new location, relative to some file or directory.
     *
     * @param rLeafName  the leaf name of the file/dir to find
     * @param rParentOrSibling  where to look for it
     */
    virtual void SetPath(const std::string& rLeafName, const FileFinder& rParentOrSibling);

    /**
     * @return true if this FileFinder has been given a path.
     */
    bool IsPathSet() const;

    /**
     * @return true if we exist (as either a file or a directory).
     */
    bool Exists() const;

    /**
     * @return true if we are pointing at a file
     */
    bool IsFile() const;

    /**
     * @return true if we are pointing at a directory
     */
    bool IsDir() const;

    /**
     * @return true if this is a file of size zero or
     * if this is a folder, whether it contains no non-hidden items.
     * If this doesn't exist, throws.
     */
    bool IsEmpty() const;

    /**
     * @return the absolute path to this file/dir.
     *
     * If this is a directory that exists (at the instant of this call), the absolute path is
     * guaranteed to end in a '/'.  Otherwise, the path is guaranteed not to end in a '/'.
     */
    std::string GetAbsolutePath() const;

    /**
     * @return true if this file/dir is newer than another file/dir.
     * Compares modification times.
     *
     * @param rOtherEntity  the entity to test against.
     */
    bool IsNewerThan(const FileFinder& rOtherEntity) const;

    /**
     * @return the leaf name of this file or directory.
     *
     * i.e. the individual file or directory name and none of the preceding folders on its path.
     */
    std::string GetLeafName() const;

    /**
     * @return the leaf name of this file or directory, with any file extension removed.
     *
     * i.e. the individual file or directory name and none of the preceding folders on its path.
     */
    std::string GetLeafNameNoExtension() const;

    /**
     * @return the extension of the leaf name of this file or directory, if any.
     * The '.' will be included in the extension if an extension exists.
     */
    std::string GetExtension() const;

    /**
     * @return a finder for the folder containing this file or directory.
     */
    FileFinder GetParent() const;

    /**
     * @return the relative path to this finder from another.  Throws if this is not found under rBasePath.
     *
     * @param rBasePath  where the returned path should be relative to
     */
    std::string GetRelativePath(const FileFinder& rBasePath) const;

    /**
     * Copy this file or folder (recursively in the latter case) to the given destination.
     * If the destination is a folder that exists, the source will be copied with the same
     * name inside that folder.  Otherwise the source will be copied with the given destination
     * name.
     *
     * If the source is a file and the destination is a file that exists it will be removed prior to copying.
     * If the source is a folder and the destination is a file that exists then an error is thrown.
     *
     * @param rDest  where to copy to
     * @return  a finder for the copied entity
     */
    FileFinder CopyTo(const FileFinder& rDest) const;

    /**
     * Recursively remove this file or folder.
     * Since this is a potentially very dangerous operation, only locations under the Chaste
     * test output folder may be removed.
     *
     * Only folders created by an OutputFileHandler, or the contents of such a folder, may be
     * deleted (folders that have .chaste_deletable_folder present).
     *
     * If you need to delete a file or folder without .chaste_deletable_folder, then you have to use
     * DangerousRemove().
     *
     */
    void Remove() const;

    /**
     * This method will allow you to remove any file from under either
     *  * CHASTE_TEST_OUTPUT
     *  or
     *  * the source tree
     * (but not elsewhere).
     *
     * For this reason it is a very dangerous operation and should not be used if Remove could be instead.
     *
     * BEWARE: if you have managed to set CHASTE_TEST_OUTPUT to "/" this could wipe your system!
     */
    void DangerousRemove() const;

    /**
     * @return a list of files in this folder matching a simple glob pattern.
     * This method must be called on a FileFinder that points at a folder, and the pattern
     * will be matched against file (or folder) names in that folder.  The pattern can use
     * a subset of shell-style glob syntax.  A '?' anywhere in the string matches any single
     * character at that position.  A '*' may be used at the start or end of the string to
     * match any number of leading or trailing characters, respectively.
     * Hidden files (names starting with a '.') will never be matched. Returns a sorted alphabetical list.
     *
     * @param rPattern  the pattern to match names against
     */
    std::vector<FileFinder> FindMatches(const std::string& rPattern) const;

    /**
     * @return true if the path is absolute.
     *
     * @param rPath  the path to test
     */
    static bool IsAbsolutePath(const std::string& rPath);

    /**
     * Replace any spaces in a path or filename with underscores.
     *
     * @param rPath  a path or file name
     */
    static void ReplaceSpacesWithUnderscores(std::string& rPath);

    /**
     * Replace any underscores in a path or filename with spaces (for making titles etc.).
     *
     * @param rPath  a path or file name
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

    /**
     * Provide a sort operator to get a logical ordering from FindMatches
     * it orders by alphabetical (or ASCII really).
     * @param otherFinder  Another FileFinder
     * @return Whether this FileFinder is earlier in the alphabetical ordering than otherFinder
     * */
    bool operator<(const FileFinder& otherFinder) const;

private:
    /** The absolute path to our file. */
    std::string mAbsPath;

    /** Whether to fake one of the fixed paths, e.g. ChasteSourceRoot. */
    static bool msFaking;

    /** Which path to fake. */
    static RelativeTo::Value msFakeWhat;

    /** The fake value of the faked path. */
    static std::string msFakePath;

    /**
     * This is code common to Remove() and DangerousRemove(). Should remain private and not to be called from elsewhere.
     * Remove() is only allowed to delete things with a .chaste_deletable_folder in the testoutput directory.
     *
     * DangerousRemove() is allowed to delete anything in the source or testoutput directories.
     *
     * @param dangerous  whether we are doing a dangerous remove.
     */
    void PrivateRemove(bool dangerous = false) const;
};

#endif /*FILEFINDER_HPP_*/
