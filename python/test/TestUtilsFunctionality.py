
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

import os
import unittest
import sys
import filecmp

class TestUtilsFunctionality(unittest.TestCase):
    
    def TestAddingVtuAnnotationsMonodomain2d(self):
        #without suffices (input and output file names)
        original_vtu1 = 'python/test/data/input/Monodomain2d_before_annotations'
        output_vtu1 = os.path.join(CHASTE_TEST_OUTPUT, 'MonodomainAddingAnnotationsTest')
        # Exceptional behaviour - do not wipe out the original file
        rc = os.system('python/utils/AddVtuTimeAnnotations.py ' + original_vtu1 + ' ' + original_vtu1)
        self.assertNotEqual(rc, 0) # No success: Error happened
        rc = os.system('python/utils/AddVtuTimeAnnotations.py ' + original_vtu1 + ' ' + output_vtu1)
        self.assertEqual(rc, 0)
        self.assertTrue(filecmp.cmp(output_vtu1+'.vtu','python/test/data/output/Monodomain2d_after_annotations.vtu'))
        
        #with suffices
        original_vtu2 = 'python/test/data/input/Monodomain2d_before_annotations.vtu'
        output_vtu2 = os.path.join(CHASTE_TEST_OUTPUT, 'MonodomainAddingAnnotationsTest2.vtu')
        rc = os.system('python/utils/AddVtuTimeAnnotations.py ' + original_vtu2 + ' ' + output_vtu2)
        self.assertEqual(rc, 0)
        self.assertTrue(filecmp.cmp(output_vtu2,'python/test/data/output/Monodomain2d_after_annotations.vtu'))

    def TestAddingVtuAnnotationsBidomain(self):
        original_vtu = 'python/test/data/input/Bidomain3d_before_annotations.vtu'
        output_vtu = os.path.join(CHASTE_TEST_OUTPUT, 'BidomainAddingAnnotationsTest.vtu')
        rc = os.system('python/utils/AddVtuTimeAnnotations.py ' + original_vtu + ' ' + output_vtu)
        self.assertEqual(rc, 0)
        self.assertTrue(filecmp.cmp(output_vtu,'python/test/data/output/Bidomain3d_after_annotations.vtu'))
        
    def TestAddingVtuAnnotationsAirway(self):
        original_vtu = 'python/test/data/input/ThreeBifurcations_before_annotations.vtu'
        output_vtu = os.path.join(CHASTE_TEST_OUTPUT, 'ThreeBifurcationsAddingAnnotationsTest.vtu')
        rc = os.system('python/utils/AddVtuTimeAnnotations.py ' + original_vtu + ' ' + output_vtu)
        self.assertEqual(rc, 0)
        self.assertTrue(filecmp.cmp(output_vtu, 'python/test/data/output/ThreeBifurcations_after_annotations.vtu'))
        
    def TestAddingVtuAnnotationsParallelPieces(self):
        original_pvtu = 'python/test/data/input/monodomain3d.pvtu'
        output_base_name = os.path.join(CHASTE_TEST_OUTPUT, 'MonodomainAddingAnnotations')
        output_pvtu = output_base_name+'.pvtu'
        rc = os.system('python/utils/AddVtuTimeAnnotations.py ' + original_pvtu + ' ' + output_pvtu)
        self.assertEqual(rc, 0)
        self.assertTrue(filecmp.cmp(output_pvtu,'python/test/data/output/MonodomainAddingAnnotations.pvtu'))
        self.assertTrue(filecmp.cmp(output_base_name+'_0.vtu','python/test/data/output/MonodomainAddingAnnotations_0.vtu'))
        self.assertTrue(filecmp.cmp(output_base_name+'_1.vtu','python/test/data/output/MonodomainAddingAnnotations_1.vtu'))
        
    def TestAddingVtuAnnotationsBidomainWithApdMap_Fails(self):
        #The script returns a file that seg faults in Paraview
        original_vtu = 'python/test/data/input/bidomain2d_with_apd_map.vtu'
        output_vtu = os.path.join(CHASTE_TEST_OUTPUT, 'bidomain2d_with_apd_map_output.vtu')
        rc = os.system('python/utils/AddVtuTimeAnnotations.py ' + original_vtu + ' ' + output_vtu)
        self.assertEqual(rc, 0)
        self.assertTrue(filecmp.cmp(output_vtu,'python/test/data/output/bidomain2d_with_apd_map_output.vtu'))
        #Extra file from post-processing_map_output_Apd_90_minus_30_Map.
        output_post_vtu = os.path.join(CHASTE_TEST_OUTPUT, 'bidomain2d_with_apd_map_output_Apd_90_minus_30_Map.vtu')
        self.assertTrue(filecmp.cmp(output_post_vtu,'python/test/data/output/bidomain2d_with_apd_map_output_Apd_90_minus_30_Map.vtu'))


    def TestConvertBinaryAttributesForChaste_3_2(self):
        # Original (pre version 3.2) has unsigned attributes (in binary and ascii files)
        # Newer Chaste has double attributes and can't cope with old binary file (each data item is 4 bytes shorter than expected)

        rc = os.system('python/utils/ConvertBinaryElementAttributes.py')
        self.assertEqual(rc, 256) #Not enough arguments - return 1

        rc = os.system('python/utils/ConvertBinaryElementAttributes.py SomeFile.node output.node')
        self.assertEqual(rc, 256*2) #Not correct extension - return 2

        rc = os.system('python/utils/ConvertBinaryElementAttributes.py mesh/test/data/simple_cube.ele output.ele')
        self.assertEqual(rc, 256*3) #An ascii mesh file - return 3

        converted = 'python/test/data/output/simple_cube_binary.ele'
        rc = os.system('python/utils/ConvertBinaryElementAttributes.py ' + converted +' output.ele')
        self.assertEqual(rc, 256*4) #A file which is already in the new format - return 4

        original = 'python/test/data/input/simple_cube_binary.ele'
        output = os.path.join(CHASTE_TEST_OUTPUT, 'simple_cube_binary.ele')

        rc = os.system('python/utils/ConvertBinaryElementAttributes.py ' +  original + ' ' +  output)
        # Byte-for-byte comparion
        self.assertTrue(filecmp.cmp(output, converted))

                     
        
