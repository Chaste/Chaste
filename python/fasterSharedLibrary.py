
"""Copyright (c) 2005-2012, University of Oxford.
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
