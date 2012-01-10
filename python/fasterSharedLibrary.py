
"""Copyright (C) University of Oxford, 2005-2012

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
"""


####################################################
# override SharedLibrary to only update the libs in '#shlib/' if the symbols
# change.  This is to avoid rebuilding binaries when a shared library has
# changes.
# This would be MUCH simpler to implement if we could override the
# signature to used for the SharedLibrary node itself directly. That is
# certainly possible, but would rely on the internal structure of SCons.

def fasterSharedLibrary(env, library, sources, **args):
    # SCons version compatibility
    if type(library) != type([]):
        library = [library]
    # use the 'quicker' shallow copy method!
    envContentSig = env.Clone()
    envContentSig.TargetSignatures('content')

    cat = env.OriginalSharedLibrary(library, sources)

    # copy all the latest libraries to ONE directory..
    # for our convenience. Could modify the above to
    # build directly to this dir instead.
    catLib = env.Install('#lib', cat) #the CURRENT lib dir

    # now generate the 'interface' file, using the
    # content signature for its target
    catIF = envContentSig.Command(
        '%s.if' % library[0],
        catLib,
        'nm --extern-only $SOURCES | cut -c 12- | sort > $TARGET')

    # install command to copy lib to shlib, where the link
    # actually occurs.  Explicitly make this depend only on
    # the IF file, which has a target content signature.
    # ie only if the Global Symbol list changes, is copied and this the
    # Programs it relinked.
    catLink = env.Command(
        '#linklib/${SHLIBPREFIX}%s${SHLIBSUFFIX}' % library[0],
        '',
        Copy('$TARGET', str(catLib[0])))
    env['CHASTE_LIBRARIES'][library[0]] = catLink[0]

    envContentSig.Depends(catLink, catIF)

    return cat

# declaring OriginalSharedLibrary is a bit marginal.  Probably should use
# a functor style object so we can store it in side the object?
#env['BUILDERS']['OriginalSharedLibrary'] = env['BUILDERS']['SharedLibrary']
#env['BUILDERS']['SharedLibrary'] = fasterSharedLibrary
