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

#include "SteadyStateExtrapolationRunner.hpp"

void SteadyStateExtrapolationRunner::RunToSteadyStateImplementation(){
  mpModel->SetMinimalReset(true);
  for(mNumEvaluations  = 0; mNumEvaluations < mMaxNumPaces; mNumEvaluations++){
    if(RunPace() == true){
      //We have reached the threshold
      return;
    }
  }
  mpModel->SetMinimalReset(false);
  return;
}


bool SteadyStateExtrapolationRunner::RunPace(){
    bool extrapolated = false;

    extrapolated = ExtrapolateStates();

    if(!extrapolated){
    /*Solve in two parts*/
      try{
	mpModel->SolveAndUpdateState(0, mDuration);
	mpModel->SolveAndUpdateState(mDuration, mPeriod);
      }
      catch(Exception &e){
	std::cout << "RunPace() failed - returning to old state_variables\n";
	mStateVariables = mSafeStateVariables;
	mpModel->SetStateVariables(mStateVariables);
	mMrmsBuffer.clear();
	mStatesBuffer.clear();
	return false;
      }
      std::vector<double> new_state_variables = mpModel->GetStdVecStateVariables();
      mStatesBuffer.push_back(new_state_variables);
      double current_mrms = CalculateMrms(new_state_variables, mStateVariables);
      mMrmsBuffer.push_back(current_mrms);
      mStateVariables = new_state_variables;
      if(current_mrms < mThreshold){
	return true;
      }
    }
    else{
      mStatesBuffer.clear();
      mMrmsBuffer.clear();
    }
    return false;
  }
     
double SteadyStateExtrapolationRunner::CalculateMrms(std::vector<double> A, std::vector<double> B){
    double norm = 0;
    
    for(unsigned int i=0; i < A.size(); i++){
      double a = A[i];
      double b = B[i];
      norm += pow((a - b)/(1 + abs(a)), 2);   
    }
    return sqrt(norm/A.size());
  }

bool SteadyStateExtrapolationRunner::ExtrapolateState(unsigned int state_index){
 
    /* Calculate the log absolute differences of the state and store these in y_vals. Store the corresponding x values in x_vals*/
    std::vector<double> y_vals, x_vals;  
    std::vector<double> state;
    state.reserve(mStatesBuffer.size());

    for(unsigned int i = 0; i < mBufferSize; i++){
      state.push_back(mStatesBuffer[i][state_index]);
    }
    
    for(unsigned int i = 0; i < mBufferSize - 1; i++){
      double tmp = abs(state[i] - state[i+1]);
      if(tmp != 0){
	y_vals.push_back(log(tmp));
	x_vals.push_back(i);
      }
    }
    const double pmcc = CalculatePMCC(x_vals, y_vals);

    if(pmcc > -0.8){
      return false;
    }
    
    /*Compute the required sums*/
    
    double sum_x = 0, sum_y = 0, sum_x2 = 0, sum_xy = 0;
    const unsigned int N = x_vals.size();

    if(N<=2){
      return false;
    }
    
    for(unsigned int i = 0; i < N; i++){
      sum_x += x_vals[i];
      sum_y += y_vals[i];
      sum_x2+= x_vals[i]*x_vals[i];
      sum_xy+= x_vals[i]*y_vals[i];
    }

    /*alpha and beta are the simple linear regression parameters*/
    
    const double beta  = (N*sum_xy - sum_x*sum_y) / (N*sum_x2 - sum_x*sum_x);

    const double alpha = (sum_y*sum_x2 - sum_x*sum_xy) / (N*sum_x2 - sum_x*sum_x);

    if(pmcc>-0.75){  //return the unchanged value if there is no negative correlation (PMCC > -0.0.75 or PMCC = NAN) 
      return false;
    }
       
    if(exp(alpha) < 100*TolRel*state.back()){
      //The difference will be about as small as solver tolerances so there is no point going any further
      return false;
    }
    
    if(beta > 0){
      //The difference is increasing, this contradicts the assumption that we have exponential decay.
      return false;
    }
    

    const double tau = -1/beta;
    double change_in_variable =  abs(mExtrapolationCoefficient * exp(alpha - mBufferSize/tau + 1/tau) / (exp(1/tau) - 1));

    /*Is V(t) increasing or decreasing?*/    
    
    if(state.back() - state.front() < 0)
      change_in_variable = - change_in_variable;
    
    double new_value = state.back() + change_in_variable; 
    
    if(std::isfinite(new_value)){
      mStateVariables[state_index] = new_value;
      return true;
    }
    else{
      return false;
    }   
}
  
bool SteadyStateExtrapolationRunner::ExtrapolateStates(){
    if(mJumps>=mMaxJumps)
      return false;
    if(!mMrmsBuffer.full())
      return false;
    double mrms_pmcc = CalculatePMCC(mMrmsBuffer);
    std::vector<double> new_state_variables;
    bool extrapolated = false;
    std::string model_name = mpModel->GetSystemInformation()->GetSystemName();
    if(mrms_pmcc < -0.95){
      mSafeStateVariables = mStateVariables;
      
      for(unsigned int i = 0; i < mStateVariables.size(); i++){
       if(ExtrapolateState(i))
	  extrapolated = true;
      }
      
      if(extrapolated){
	mMrmsBuffer.clear();
	mStatesBuffer.clear();
	mpModel->SetStateVariables(mStateVariables);
	mJumps++;
      }
    }
    return extrapolated;
  }
  
