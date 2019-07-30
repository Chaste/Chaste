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

#ifndef DEBUG_HPP_
#define DEBUG_HPP_

#include <iostream>
#include <cassert>
#include <sstream>
#include <string>
#include "PetscTools.hpp"

/**
 * @file
 * A bunch of useful macros for debugging.  These all output information to stdout as
 * lines starting "DEBUG: ".  When running in parallel the process number is also
 * displayed.  Each line is flushed as soon as it is written.
 *
 * @note Use of these should be removed from source code when committing.
 */

/** Can be set as approximate memory footprint in Mb at marker using the macro MARK_MEMORY */
extern double eMemoryAtMarker;

/** @return a 'header' which can be printed in each debug output line */
std::string FormDebugHead();

/**
 * In extremis, print enough of the stack to see how we got to this point (e.g. at EXCEPTION time).
 * This is used with the macro STACK
 */
void PrintTheStack();

/**
 * Allows use of PRINT_MEMORY to give results relative to the footprint at this call.
 *
 * Set sMemoryAtMarker to the current approximate memory footprint.  This value is subtracted
 * from the memory footprint if the macro PRINT_MEMORY is called.
 *
 * This is used with the macro MARK_MEMORY
 */
void MarkMemory();

/**
 * Set sMemoryAtMarker to zero meaning subsequent use of PRINT_MEMORY will give an
 * absolute memory footprint, rather than relative to a fixed value.
 *
 * This is used with the macro UNMARK_MEMORY
 */
void UnmarkMemory();

/**
 * Print the approximate memory footprint.  This is the current memory footprint minus
 * sMemoryAtMarker, which will be zero if the MARK_MEMORY macro has not been used.
 *
 * This is used with the macro PRINT_MEMORY
 */
void PrintMemory();

/**
 * Print a debug message.
 * @param stuff  what to print (can be a variable, or e.g. a << " " << b)
 */
#define TRACE(stuff) std::cout << FormDebugHead() << stuff << std::endl << std::flush;

/** Print some trace containing the file name and line number. */
#define MARK std::cout << FormDebugHead() <<  __FILE__ << " at line " << __LINE__ << std::endl << std::flush;
/** Print some trace containing the file name and line number, but only on the matching process*/
#define MARK_ON_PROCESS(proc) if (PetscTools::GetMyRank()==proc) std::cout << FormDebugHead() <<  __FILE__ << " at line " << __LINE__ << std::endl << std::flush;
/** Print the name and value of the given variable.
 * @param var */
#define PRINT_VARIABLE(var) std::cout << FormDebugHead() << #var " = " << var << std::endl << std::flush;

/** Print the name and value of the given variables.
 * @param var1
 * @param var2
 */
#define PRINT_2_VARIABLES(var1,var2) std::cout << FormDebugHead() << #var1 " = " << var1 << ", " \
    #var2 " = " << var2 << std::endl << std::flush;

/** Print the name and value of the given variables.
 * @param var1
 * @param var2
 * @param var3
 */
#define PRINT_3_VARIABLES(var1,var2,var3) std::cout << FormDebugHead() << #var1 " = " << var1 << ", " \
    #var2 " = " << var2 << ", " #var3 " = " << var3 << std::endl << std::flush;

/** Print the name and value of the given variables.
 * @param var1
 * @param var2
 * @param var3
 * @param var4
 */
#define PRINT_4_VARIABLES(var1,var2,var3,var4) std::cout << FormDebugHead() << #var1 " = " << var1 << ", " \
    #var2 " = " << var2 << ", " #var3 " = " << var3 << ", " \
    #var4 " = " << var4 << std::endl << std::flush;

/** Print the name and value of the given variables.
 * @param var1
 * @param var2
 * @param var3
 * @param var4
 * @param var5
 */
#define PRINT_5_VARIABLES(var1,var2,var3,var4,var5) std::cout << FormDebugHead() << #var1 " = " << var1 << ", " \
    #var2 " = " << var2 << ", " #var3 " = " << var3 << ", " \
    #var4 " = " << var4 << ", " #var5 " = " << var5 <<std::endl << std::flush;

/** Quit (assert(0)) on the n-th time this line is reached, for the given n.
 * @param n */
#define QUIT_AFTER_N_VISITS(n) { static unsigned counter=0; if (++counter==(n)) {TRACE("User-forced quit."); assert(0);} }

/** Print how many times this line has been reached, everytime it is reached.
 * @param message  message to include in brackets */
#define HOW_MANY_TIMES_HERE(message) { \
    static unsigned counter=1; \
    std::cout << FormDebugHead() << "Num times here (" << message << "): " << counter++ << std::endl << std::flush; }

/** Prints the given message, but only from the n-th time that line is reached, for the given n.
 * @param stuff  what to print
 * @param n */
#define TRACE_FROM_NTH_VISIT(stuff,n) { \
    static unsigned counter=0; \
    if (++counter>=(n)) {TRACE(stuff<<" (visit "<<counter<<")");} }

/** Display a std::vector.
 * @param v */
#define PRINT_VECTOR(v) \
    { std::cout << FormDebugHead() << #v " = {"; \
      for (unsigned _i=0; _i<v.size(); _i++) { \
          std::cout << (_i==0?"":",") << v[_i]; } \
      std::cout << "}" << std::endl << std::flush; }

/** (For debugging on a seriously large number of processes)
 * Show the mark that we have reached this line of code, but do it in process order, which
 * is important if we suspect that one or more processes are missing.
 *
 * This will also isolate all the MARKs which appear for a particular line of code, since all process have to synchronise at
 * the beginning of the round-robin.
 *
 * Note that this will change the parallel behaviour of the code -- if one process really is missing then the first barrier inside
 * the round-robin will cause deadlock.
 */
#define MARK_IN_ORDER PetscTools::BeginRoundRobin(); MARK; PetscTools::EndRoundRobin();

/**
 * In extremis, print enough of the stack to see how we got to this point (e.g. at EXCEPTION time)
 */
#define STACK PrintTheStack();

/**
 * Allows use of PRINT_MEMORY to give results relative to the footprint at this call.
 *
 * Set sMemoryAtMarker to the current approximate memory footprint.  This value is subtracted
 * from the memory footprint if the macro PRINT_MEMORY is called.
 */
#define MARK_MEMORY MarkMemory();

/**
 * Set sMemoryAtMarker to zero meaning subsequent use of PRINT_MEMORY will give an
 * absolute memory footprint, rather than relative to a fixed value.
 */
#define UNMARK_MEMORY UnmarkMemory();

/**
 * Print the approximate memory footprint.  This is the current memory footprint minus
 * sMemoryAtMarker, which will be zero if the MARK_MEMORY macro has not been used.
 */
#define PRINT_MEMORY PrintMemory();

#endif /*DEBUG_HPP_*/
