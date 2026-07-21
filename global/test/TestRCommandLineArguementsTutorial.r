#Copyright (c) 2005-2026, University of Oxford.
#All rights reserved.

#University of Oxford means the Chancellor, Masters and Scholars of the
#University of Oxford, having an administrative office at Wellington
#Square, Oxford OX1 2JD, UK.
#
#This file is part of Chaste.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#* Redistributions of source code must retain the above copyright notice,
#  this list of conditions and the following disclaimer.
#* Redistributions in binary form must reproduce the above copyright notice,
#  this list of conditions and the following disclaimer in the documentation
#  and/or other materials provided with the distribution.
#* Neither the name of the University of Oxford nor the names of its
#  contributors may be used to endorse or promote products derived from this
#  software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
#LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
#GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Note that this tutorial will require both testthat and tidyverse to be installed
# For help installing and managing R libraries the offical guide can be found here : https://cran.r-project.org/doc/manuals/r-patched/R-admin.html
library(testthat)

# Note : Before executing this tutorial or any other instance of Chaste tests we first need to build Chaste
# Instructions on on your first build of Chast can be found here : https://chaste.github.io/docs/user-guides/cmake-first-run/
# Within this tutorial we will show how to execute Chaste simulations/tests from the command line using R
# You can simillarly find information on executing code from Bash and Python : https://chaste.github.io/docs/user-tutorials/commandlinearguments/
# First we must define where our test can be found within within our build directory
# Within the context of this Tutorial we assume that Chaste has been built within the User's home directory
# This path will need to be changed if your build folder for Chaste is not within the home directory
path_to_chaste_build <- file.path("/home", Sys.getenv("USER"),"/build/global/test/TestCommandLineArgumentsTutorial")

# The underlying Chaste test we are attempting to execute can be found here : https://github.com/Chaste/Chaste/blob/develop/global/test/TestCommandLineArgumentsTutorial.hpp

# For our first example we will pass command line parameters to Chaste using the test_that library with a small simple loop
# The results of each loop will then be output to terminal
test_that("test_command_line_default_tutorial", {
  N <- 2
  L <- 3
  M <- 4
  
  # R's seq(start, end) is inclusive of both ends
  for (i in 0:N) {
    for (j in 1:L) {
      for (k in 2:M) {
        # system2 is the R equivalent to subprocess.run which is found within the Python equivalent tutorial
        # args is a character vector of the flags and values , these are what will be passed to our Chaste test
        args <- c("-opt1", as.character(i), 
                  "-opt2", as.character(j), 
                  "-opt3", as.character(k))
        
        status <- system2(path_to_chaste_build, args = args)
        
        # check=True equivalent: stop if the command fails
        expect_equal(status, 0, info = paste("Execution failed for i,j,k:", i, j, k))
      }
    }
  }
})

# For our second example we will pass a vector of values at command line to Chaste using the test_that library with a small simple loop
# The vector components will then be summed and output to terminal

test_that("test_command_line_double_tutorial", {
  N <- 2000
  L <- 3000
  M <- 4000
  
  # Python's range(start, stop, step) -> R's seq(from, to, by)
  for (i in seq(1, N, by = 1000)) {
    idouble <- i / 120.0
    for (j in seq(1001, L, by = 1000)) {
      jdouble <- j / 150.0
      for (k in seq(2001, M, by = 1000)) {
        kdouble <- k / 180.0
        
        args <- c("--my-vector-of-arguments", 
                  as.character(idouble), 
                  as.character(jdouble), 
                  as.character(kdouble))
        
        status <- system2(path_to_chaste_build, args = args)
        expect_equal(status, 0)
      }
    }
  }
})
