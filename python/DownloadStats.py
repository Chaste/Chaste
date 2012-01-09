#!/usr/bin/env python

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
