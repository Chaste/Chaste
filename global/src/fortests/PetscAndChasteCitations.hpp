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


/*
 * If using older an older PETSc (pre-3.2) include some citations here from a more recent version.
 * (The following are from 3.5.2)
 */
static PetscBool PetscCite1 = PETSC_FALSE;
const char PetscCitation1[] = "@TechReport{petsc-user-ref,\n"
                              "  Author = {Satish Balay and Shrirang Abhyankar and Mark F. Adams and Jed Brown and Peter Brune\n"
                              "            and Kris Buschelman and Victor Eijkhout and William D. Gropp\n"
                              "            and Dinesh Kaushik and Matthew G. Knepley\n"
                              "            and Lois Curfman McInnes and Karl Rupp and Barry F. Smith\n"
                              "            and Hong Zhang},\n"
                              "  Title = {{PETS}c Users Manual},\n"
                              "  Number = {ANL-95/11 - Revision 3.5},\n"
                              "  Institution = {Argonne National Laboratory},\n"
                              "  Year = {2014}\n"
                              "}\n";
static PetscBool PetscCite2 = PETSC_FALSE;
const char PetscCitation2[] = "@InProceedings{petsc-efficient,\n"
                              "  Author = {Satish Balay and William D. Gropp and Lois Curfman McInnes and Barry F. Smith},\n"
                              "  Title = {Efficient Management of Parallelism in Object Oriented Numerical Software Libraries},\n"
                              "  Booktitle = {Modern Software Tools in Scientific Computing},\n"
                              "  Editor = {E. Arge and A. M. Bruaset and H. P. Langtangen},\n"
                              "  Pages = {163--202},\n"
                              "  Publisher = {Birkh{\\\"{a}}user Press},\n"
                              "  Year = {1997}\n"
                              "}\n";


/* Main Chaste citation */
static PetscBool ChasteCite = PETSC_FALSE;
const char ChasteCitation[] = "@article{mirams2013chaste,\n"
                              "  author    = {Mirams, G.R. and Arthurs, C.J. and Bernabeu, M.O. and Bordas, R. and Cooper, "
                              "J. and Corrias, A. and Davit, Y. and Dunn, S-J. and Fletcher, A.G. and Harvey, D.G. and "
                              "Marsh, M.E. and Osborne, J.M. and Pathmanathan, P. and Pitt-Francis, J. and Southern, J. "
                              "and Zemzemi, N. and Gavaghan, D.J.},\n"
                              "  title     = {Chaste: an open source C++ library for computational physiology and biology},\n"
                              "  journal   = {PLoS computational biology},\n"
                              "  volume    = {9},\n"
                              "  number    = {3},\n"
                              "  pages     = {e1002970},\n"
                              "  year      = {2013},\n"
                              "  publisher = {Public Library of Science}\n"
                              "}\n";