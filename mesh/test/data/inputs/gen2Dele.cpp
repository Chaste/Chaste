/*

Copyright (C) University of Oxford, 2005-2011

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
