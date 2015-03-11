// File:  getDiscoverySignificance.cc
// Glen Cowan
// RHUL Physics

// Returns significance Z with which hypothesized value of mu=0 is rejected; 
// Z = Phi^-1(1 - p), p = p-value of mu, Phi^-1 = quantile of standard Gaussian. 
// 
// Inputs:     
// n     number of events seen in data, e.g. generate from Poisson(mu*s+b)
// s     expected number of signal events
// m     vector of numbers of events seen in subsidiary measurements.
// tau   defined by m_i ~ Poisson(tau_i*b_i)

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cmath>

#include "SigCalc.h"

double getDiscoverySignificance(double n, double s,
       vector<double> m, vector<double> tau){

  SigCalc* sc = new SigCalc (n, s, m, tau, 0);
  double q0Val = sc->q0();
  double Z = sqrt(q0Val);

  delete sc;
  return Z;

}
