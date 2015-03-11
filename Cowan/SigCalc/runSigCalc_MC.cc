// File:  runSigCalc.cc
// Glen Cowan
// RHUL Physics
// 14 April 2008

// Simple driver program for the function getSignificance, which 
// returns the significance with which one can reject a hypothesized
// signal strength parameter mu.  (For discovery, one rejects mu = 0.)

// Input:  runSigCalc reads in a text file with the input data.  
// Lines beginning with # are comments.
// The first (non-comment) line contains the target luminosity.
// The second (non-comment) line contains the accepted number of
// signal events from MC and the luminosity of the signal MC sample.
//
// On subsequent lines the user lists for each background component the
// accepted number of background events m and the effective luminosity 
// L_MC used to generate the background sample.

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <TMath.h>
#include <TRandom3.h>
#include <TH1D.h>
#include <TFile.h>
#include "SigCalc.h"
#include "getDiscoverySignificance.h"
#include "getExclusionSignificance.h"

using namespace std;

int main(int argc, char **argv) {

  TFile* histFile = new TFile("SigCalc.root", "recreate");
  TH1D* h_qmu = new TH1D("h_qmu", "qmu",  200,  0, 40.);

// Read input file

  string inputFileName;
  if ( argc == 1 ) {
    cout << "Enter name of input file: ";
    cin >> inputFileName;
    cin.ignore();         // clear buffer
  }
  else {
    inputFileName = argv[1];
  }

  ifstream f;
  f.open( inputFileName.c_str() );
  if ( f.fail() ){
    cout << "Sorry, couldn't open input file" << endl;
    exit(1);
  }

  double nObs;
  double s;
  vector<double> mObs;
  vector<double> tau;
  istringstream instream;       // Declare an input string stream
  string line;
  int lineNum = 0;
  while ( getline(f, line) ) {
    bool useLine = !( line.substr(0,1) == "#" || line.substr(0,1) == "!");
    if ( useLine ) {
      instream.clear();          // Reset from possible previous errors
      instream.str(line);        // Use line as source of input
      lineNum++;      
      if ( lineNum == 1 ) {
        instream >> nObs;           // number of observed events
      }
      else if ( lineNum == 2 ) {
        instream >> s;    // expected number of events in nominal signal model
      }
      else {
        double m_i, tau_i;
        instream >> m_i >> tau_i;      
        mObs.push_back(m_i);
        tau.push_back(tau_i);
      }
    }
  }
  int numBkg = mObs.size();

  cout << "nObs = " << nObs << endl;
  cout << "s = " << s << endl;
  double btotHat = 0.;
  for (int i=0; i<mObs.size(); i++){
    btotHat += mObs[i]/tau[i];
    cout << "m[" << i << "]  = " << mObs[i] << ",   tau[" << i << "]  = " 
         << tau[i] << endl;
  }
  cout << "observed estimated total background = " << btotHat << endl;
  cout << endl;

  int seed = 12345;
  TRandom3* ran = new TRandom3(seed);

  // Loop over values of mu, compute p-value using asymptotic formula and MC

  cout << "mu        pmu (asymptotic)       p_mu (MC)" << endl;

  double muTestMin = 0.1;
  double muTestMax = 2.0;
  int numMuVal = 20;
  for (int i=0; i<numMuVal; i++){

    double muTest = (muTestMax - muTestMin)*static_cast<double>(i) / 
                    static_cast<double>(numMuVal - 1)  +  muTestMin;

    SigCalc* sc = new SigCalc(nObs, s, mObs, tau, 0);
    double muHatObs = (nObs - btotHat) / s;       // start value in fit
    vector<double> bHatObs(numBkg);
    vector<double> bHatHatObs(numBkg);
    double qmuObs = sc->qmu(muTest, muHatObs, bHatObs, bHatHatObs);
    double ZExclObs = sqrt(qmuObs);
    double pmu = 1. - TMath::Freq(ZExclObs);
    delete sc;

    double bHatHatObsTot = 0.;
    for (int j=0; j<numBkg; j++){
      bHatHatObsTot += bHatHatObs[j];
    }

    // Now generate events using bHatHatObs values 

    int numEvt = 0;
    int numRej = 0;
    int numRejMin = 100;
    bool generateEvents = true;
    while ( generateEvents ) {

      numEvt++;
      int n = ran->Poisson(muTest*s + bHatHatObsTot);
      double nd = static_cast<double>(n);
      vector<double> md;
      md.clear();
      for (int k=0; k<numBkg; k++){
        int m = ran->Poisson(tau[k]*bHatHatObs[k]);
        md.push_back(static_cast<double>(m));
      }
      SigCalc* sc = new SigCalc(nd, s, md, tau, 0);
      double muHat = (nd - bHatHatObsTot) / s;     // start value in fit
      vector<double> bHat(numBkg);
      vector<double> bHatHat(numBkg);
      double qmu = sc->qmu(muTest, muHat, bHat, bHatHat);
      delete sc;
      h_qmu->Fill(qmu);
      if (qmu >= qmuObs) {
        numRej++;
      }
      generateEvents = numRej < numRejMin;

    }
    double pmuMC = static_cast<double>(numRej) / 
                   static_cast<double>(numEvt);
    
    cout << muTest << "  " << pmu << "  " << pmuMC << endl;

  }


  histFile->Write();
  histFile->Close();

  return 0;

}
