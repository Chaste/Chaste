
"""Copyright (c) 2005-2024, University of Oxford.
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

import unittest
import numpy as np
import chaste
chaste.init()
import chaste.mesh


class TestChastePoint(unittest.TestCase):

    def test_construct(self):

        chaste_point = chaste.mesh.ChastePoint3(1.0, 2.0, 3.0)
        loc = chaste_point.rGetLocation()
        self.assertAlmostEqual(loc[0], 1.0, 2)
        self.assertAlmostEqual(loc[1], 2.0, 2)
        self.assertAlmostEqual(loc[2], 3.0, 2)

        point2 = chaste.mesh.ChastePoint3(np.array([4.0, 5.0, 6.0]))
        loc = point2.rGetLocation()
        self.assertAlmostEqual(loc[0], 4.0, 2)
        self.assertAlmostEqual(loc[1], 5.0, 2)
        self.assertAlmostEqual(loc[2], 6.0, 2)

        point3 = chaste.mesh.ChastePoint2(np.array([4.0, 5.0]))
        loc = point3.rGetLocation()
        self.assertAlmostEqual(loc[0], 4.0, 2)
        self.assertAlmostEqual(loc[1], 5.0, 2)


if __name__ == '__main__':
    unittest.main()