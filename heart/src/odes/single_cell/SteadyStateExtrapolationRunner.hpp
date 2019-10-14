/*

Copyright (c) 2005-2018, University of Oxford.
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


#ifndef _STEADYSTATEEXTRAPOLATIONRUNNER_HPP_
#define _STEADYSTATEEXTRAPOLATIONRUNNER_HPP_

#ifdef CHASTE_CVODE

#include "AbstractSteadyStateRunner.hpp"
#include <boost/circular_buffer.hpp>
//#include "VectorHelperFunctions.hpp"

class SteadyStateExtrapolationRunner : public AbstractSteadyStateRunner
{
  
  /**
 * This class is to get a cell model to steady state by using extrapolation.
 *
 * We use the Mixed Root Mean Square (MRMS) error to measure how far from the steady state 
 * we are. The PMCC of the MRMS errors is then used to check whether we can approxiamte the 
 * decay of the state variables with an exponential. We then use simple linear regression 
 * to find this exponential and extrapolate the variables.
 */
  
private:
  /* The number of points to use for our PMCC calculation and simple linear regression*/
  const unsigned int mBufferSize = 50;
  /*The maximum number of times we will attempt to extrapolate*/
  const unsigned int mMaxJumps = 100;
  /*The MRMS value for which we are close enough to the steady state that the absolute error in APD90 is probably less that 0.01ms*/ 
  const double mThreshold = 1.8e-07;
  /** Relative and absolute tolerances */
  const double TolRel = 1e-8;
  const double TolAbs = 1e-8;
  /** How far to extrapolate towards our approximation of the limit of the state variable. */
  double  mExtrapolationCoefficient = 0.9;
  /** The period and duration of the regular stimulus function */
  double  mPeriod;
  double  mDuration;
  /**The terminal state variables of the last pace we calculated*/ 
  std::vector<double> mStateVariables;
  /*The final state variables before the last extrapolation. This is used to reset in case an extrapolation causes an error*/
  std::vector<double> mSafeStateVariables;

  /** The stimulus used to pace the model - this must be a RegularStimulus */
  boost::shared_ptr<RegularStimulus> mpStimulus;

  /**Contains mBufferSize previous values of each state variable*/ 
  boost::circular_buffer<std::vector<double>>  mStatesBuffer;
  /**Contains mBufferSize previous values of each state vMRMS value*/ 
  boost::circular_buffer<double> mMrmsBuffer;

  /**Keeps track of the number of times we have extrapolated*/
  unsigned int mJumps = 0;

  /**Run the cell to steady state */
  virtual void RunToSteadyStateImplementation();
  
  /**
     * Calculate the MRMS error between two vectors
     *
     * @param A  The second vector
     * @param B  The first  vector
     */
  double CalculateMrms(std::vector<double> A, std::vector<double> B);

  /**
     * Calculate the PMCC of a circular_buffer of variables. This is used to calculate the PMCC of the previous MRMS values.
     *
     * @param values The values that we are calculating the PMCC of. 
     */
  
  
  double CalculatePMCC(boost::circular_buffer<double> values){
    if(values.size()<=2 || values.size() <= 2){
      return -NAN;
    }
    
    const unsigned int N = values.size();
    const double sum_x = N*(N-1)/2;
    const double sum_x2 = (N-1)*N*(2*N-1)/6;
    
  double  sum_y = 0, sum_y2 = 0, sum_xy = 0;

  for(unsigned int i = 0; i < N; i++){
    sum_y += values[i];
    sum_y2+= values[i]*values[i];
    sum_xy+= i*values[i];
  }
  double pmcc = (N*sum_xy - sum_x*sum_y)/sqrt((N*sum_x2 - sum_x*sum_x)*(N*sum_y2 - sum_y*sum_y));
  return pmcc;
  }

    /**
     * Calculate the PMCC of a circular_buffer of (x,y) pairs. This is used to calculate the PMCC of the previous MRMS values.
     *
     * @param x The x values of the (x,y) pairs
     * 
     * @param y The y values of the (x,y) pairs
 */
  
  
  double CalculatePMCC(std::vector<double> x, std::vector<double> y){
    const unsigned int N = x.size();
    if(x.size() <= 2){
      /*We don't have enough values*/
      return -NAN;
    }
    double sum_x = 0, sum_x2 = 0, sum_y = 0, sum_y2 = 0, sum_xy = 0;

    /*Calculate the required sums */
    for(unsigned int i = 0; i < N; i++){
      sum_x  += x[i];
      sum_x2 += x[i]*x[i];
      sum_y  += y[i];
      sum_y2 += y[i]*y[i];
      sum_xy += x[i]*y[i];
    }
    /*Return the PMCC*/
    return (N*sum_xy - sum_x*sum_y)/sqrt((N*sum_x2 - sum_x*sum_x)*(N*sum_y2 - sum_y*sum_y));
  }

  /** 
   * Applies the extrapolation method to one state variable if it is sensible to do so. Returns
   * true if the state variable has been extrapolated.
   *
   * @state_index The index of the state variable that we wish to extrapolate
   */
  bool ExtrapolateState(unsigned int state_index);

   /** 
   * Applies the extrapolation method to every state variable
   *
   * @state_index The index of the state variable that we wish to extrapolate
   */
  bool ExtrapolateStates();

  /** 
   * Simulates one pace and extrapolates the state variables if the PMCC
   * of the previous MRMS values is approximately -1
   */ 
  bool RunPace();
  
public:
/**
     * Constructor of a helper class for getting action potential models to steady state.
     *
     * @param pModel  The cell model to run to steady state.
     * @param period  The period for the model to be paced at.
     */
  
  SteadyStateExtrapolationRunner(boost::shared_ptr<AbstractCvodeCell> pModel, double period)
    : AbstractSteadyStateRunner(pModel), mPeriod(period){
    /*This ensures that we are using a regular stimulus*/
    mpStimulus = mpModel->UseCellMLDefaultStimulus();
    mpStimulus->SetPeriod(mPeriod);
    mDuration = mpStimulus->GetDuration();

    
    /* Set Solver Tolerances */
    mpModel->SetTolerances(TolRel, TolAbs);

    /* Set the maximum number of timesteps for CVODE to use. This prevents 
     * CV_TOO_MUCH_WORK errors */
    mpModel->SetMaxSteps(1e5);

    /* Initialise mStateVariables */
    mStateVariables = mpModel->GetStdVecStateVariables();

    /* Set the capacity of the circular buffers */
    mMrmsBuffer.set_capacity(mBufferSize);
    mStatesBuffer.set_capacity(mBufferSize);    
  };
  
}; 


#endif // CHASTE_CVODE

#endif // _STEADYSTATEEXTRAPOLATIONRUNNER_HPP_

