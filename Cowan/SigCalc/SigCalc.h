// File:  SigCalc.h
// Glen Cowan
// RHUL Physics

#ifndef SIGCALC_H
#define SIGCALC_H

#include <iostream>
#include <vector>

using namespace std;

class SigCalc {

  public: 

    SigCalc (double n, double s, vector<double> mVec, 
             vector<double> tauVec, int option);
    double n() { return m_n; }
    double s() { return m_s; }
    double m(int i) { return m_m[i]; }
    double tau(int i) { return m_tau[i]; }
    double w(int i) { return m_w[i]; }
    int numBck(){ return m_numBck; }
    double lnL(double mu, vector<double> bVec);
    double qmu(double mu, double& muHat, vector<double>& bHat,
             vector<double>& bHatHat);
    double qmu(double mu);
    double q0(double& muHat, vector<double>& bHat,
             vector<double>& bHatHat);
    double q0();

  private:

    int            m_numBck;
    double         m_n;
    double         m_s;
    vector<double> m_m;
    vector<double> m_tau;
    vector<double> m_w;

    int              m_option;

};

#endif
