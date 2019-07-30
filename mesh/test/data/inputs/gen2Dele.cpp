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


#include <iostream>

int main()
{
    int n = 0, t = 0, i = 0;

    std::cout << "200\t3\t0" << std::endl;

    while (n < 100)
    {
        if (t != 10)
        {
            // Top-left - bottom-right orientation

            std::cout << 2*n << "\t" << i << "\t" << i+1 << "\t" << i+11 << std::endl;
            std::cout << 2*n+1 << "\t" << i+1 << "\t" << i+12 << "\t" << i+11 << std::endl;

            // Bottom-left - top-right orientation

//         std::cout << 2*n << "\t" << i << "\t" << i+1 << "\t" << i+12 << std::endl;
//         std::cout << 2*n+1 << "\t" << i << "\t" << i+12 << "\t" << i+11 << std::endl;

            // Criss-cross shape - !!remember to change 200 to 400 in line above outside loop

            //         std::cout << 4*n << "\t" << i << "\t" << i+1 << "\t" << n+121 << std::endl;
            //         std::cout << 4*n+1 << "\t" << i+1 << "\t" << i+12 << "\t" << n+121 << std::endl;
            //std::cout << 4*n+2 << "\t" << i+12 << "\t" << i+11 << "\t" << n+121 << std::endl;
            //std::cout << 4*n+3 << "\t" << i+11 << "\t" << i << "\t" << n+121 << std::endl;

            t++;
            n++;
        }
        else
        {
            t = 0;
        }

        i++;
    }

    return 0;
}
