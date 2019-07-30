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

from __future__ import print_function  # print without newline character

import platform


"""
Script prints a human-readable operating system name to the terminal.

This is primarily to be executed from cmake/Modules/ChasteHostOperatingSystem.cmake, which prints the operating system
name during the CMake configure step.

"""


def main():
    platform_name = platform.platform()

    # Linux
    if platform_name.lower().startswith('linux'):
        linux_ver = platform.linux_distribution()
        return '%s %s (%s)' % (linux_ver[0], linux_ver[1], linux_ver[2])
    # E.g.                      ^ Ubuntu      ^ 14.04       ^ trusty

    # Mac
    elif platform_name.lower().startswith('darwin'):
        mac_ver = platform.mac_ver()
        mac_name = get_mac_name_from_darwin_version(mac_ver[0])
        return '%s (%s)' % (mac_name, mac_ver[0])
    # E.g.                  ^ macOS Sierra ^ 10.12.6

    # Windows
    elif platform_name.lower().startswith('windows'):
        win_ver = platform.win32_ver()
        return 'Windows %s %s' % (win_ver[0], win_ver[2])
    # E.g.                         ^ 7         ^ SP1

    # Unknown
    else:
        return 'Unknown OS'


def get_mac_name_from_darwin_version(version_str):
    if version_str.startswith('10.0.'):
        return 'OS X Cheetah'
    elif version_str.startswith('10.1.'):
        return 'OS X Puma'
    elif version_str.startswith('10.2.'):
        return 'OS X Jaguar'
    elif version_str.startswith('10.3.'):
        return 'OS X Panther'
    elif version_str.startswith('10.4.'):
        return 'OS X Tiger'
    elif version_str.startswith('10.5.'):
        return 'OS X Leopard'
    elif version_str.startswith('10.6.'):
        return 'OS X Snow Leopard'
    elif version_str.startswith('10.7.'):
        return 'OS X Lion'
    elif version_str.startswith('10.8.'):
        return 'OS X Mountain Lion'
    elif version_str.startswith('10.9.'):
        return 'OS X Mavericks'
    elif version_str.startswith('10.10.'):
        return 'OS X Yosemite'
    elif version_str.startswith('10.11.'):
        return 'OS X El Capitan'
    elif version_str.startswith('10.12.'):
        return 'macOS Sierra'
    elif version_str.startswith('10.13.'):
        return 'macOS High Sierra'
    elif version_str.startswith('10.14.'):
        return 'macOS Mojave'
    else:
        return 'Unknown Mac Version'


if __name__ == "__main__":
    print(main(), end='')
