
"""Copyright (c) 2005-2019, University of Oxford.
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

"""Host-specific configuration.

This file contains the switching on machine name which allows the core
Chaste developers to store their hostconfig information in the repository.

This whole folder should be removed when making a release.
"""


import socket
import os

def config_module():
    """Return the config module for this machine, if one exists.

    If the file 'local.py' exists, it is assumed to contain the local
    configuration.  Otherwise, we switch based on hostname to find a
    config file.  If nothing suitable is found, throw an exception.
    """

    try:
        import local
        return local
    except ImportError:
        pass # keep going below
    
    import re

    machine_fqdn = socket.getfqdn()
    name_parts = machine_fqdn.split('.', 1)
    if len(name_parts) == 1:
        name_parts.append('')
    machine_name, machine_domain = name_parts
    #print machine_fqdn, machine_name, machine_domain

    # Oxford intel compiler licence
    if machine_fqdn.endswith('ox.ac.uk'):
        os.environ['INTEL_LICENSE_FILE'] = '28519@flexlm.nsms.ox.ac.uk:' + os.environ.get('INTEL_LICENSE_FILE', '.')

    #hardy-vm is used for making the 32-bit standalone exe.
    if machine_fqdn == "hardy-vm": #and os.getuid()==1001:
        import hardyvm as conf
    # Newer ubuntu servers
    elif machine_name in ['dizzy', 'scoop', 'roley', 'travis', 'userpc58', 'muck', 'trix', 'lofty'] and machine_domain == 'cs.ox.ac.uk':
        import dizzy as conf
    #Developer machines    
    elif machine_fqdn == "heart.cs.ox.ac.uk":
        import heart as conf
    elif machine_fqdn == "beat.cs.ox.ac.uk":
        import beat as conf
    elif machine_fqdn in ["compphys0", "compphys11", "compphys00", "alberto"]:
        import alberto as conf
    elif machine_fqdn == "csu6388.cs.ox.ac.uk":
        import joe as conf
    elif machine_fqdn == "clpc287.cs.ox.ac.uk":
        import jonc as conf
    elif machine_fqdn == 'skip.cs.ox.ac.uk':
        import gary as conf
    elif machine_fqdn.startswith("caribou"):
        import caribou as conf
    elif machine_fqdn.endswith(".oerc.ox.ac.uk"):
        import oerc as conf
    elif machine_fqdn.endswith("arcus.osc.local"):
        import arcus as conf
    elif machine_fqdn.endswith("arcus.arc.local"):
        import arcus_b as conf
    elif machine_fqdn == 'clpc58.cs.ox.ac.uk':
        import kylie as conf  
    elif machine_fqdn.startswith("compute-lung"):
        import computelung as conf
    elif machine_fqdn == "tumbler.cs.ox.ac.uk":
        import tumbler as conf
    elif machine_fqdn.endswith(".cs.ox.ac.uk"):
        import comlab as conf
    elif machine_fqdn.endswith(".dtc.ox.ac.uk"):
        import dtc as conf
    elif machine_fqdn.endswith(".iceberg.shef.ac.uk"):
        import iceberg as conf
    elif machine_fqdn == "beaky-buzzard.maths.ox.ac.uk":
        import maths as conf
    elif machine_fqdn == "tweety-pie.maths.ox.ac.uk" or machine_fqdn.startswith("jochen"):
        import kursawe as conf
    elif machine_fqdn in ["tex-tucker.maths.ox.ac.uk", "fergus-laptop"]:
        import fcooper as conf
    elif machine_fqdn.startswith('alex-laptop'):
        import alexf as conf
    elif machine_fqdn.startswith('ozzy-desktop'):
        import ozzy64 as conf  
    elif re.match('eslogin[0-9]+$', machine_fqdn):
        import archer as conf
    elif machine_fqdn.startswith("hlwus011"):
        import gsk_hlwus011 as conf
    elif machine_fqdn.startswith("hlwwolx057"):
        import gsk_hlwwolx057 as conf
    else:
        raise ImportError('No config file found specifically for this machine')

    return conf
