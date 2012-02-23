#!/usr/bin/env python

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

"""
Parse the email download logs to get some information about
who has downloaded Chaste.

See also https://chaste.cs.ox.ac.uk/chaste/stats.php for
raw download numbers.
"""

email_mbox = '/users/jonc/Mail/Chaste/downloads'

# Entries to look for
keys = {'E-mail: ': 'email',
        'First name: ': 'fname',
        'Last name: ': 'lname'}
keynames = keys.values()
keynames.sort()

# Initialise
people = set()
person = {}

# Parse mailbox
for line in file(email_mbox):
    for key, var in keys.items():
        if line.startswith(key):
            person[var] = line.strip()[len(key):]
            break
    if len(person) == len(keynames):
        # Found a complete entry
        p = map(lambda k: person[k].lower(), keynames)
        ##print p
        people.add(tuple(p))
        person = {}

# Report results
people = list(people)
people.sort()
for person in people:
    print person[0], person[1], person[2]
