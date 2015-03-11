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
#include "SigCalc.h"
#include "getDiscoverySignificance.h"
#include "getExclusionSignificance.h"

using namespace std;

int main(int argc, char **argv) {

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

  double n, s;
  vector<double> m;
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
        instream >> n;           // number of observed events
      }
      else if ( lineNum == 2 ) {
        instream >> s;    // expected number of events in nominal signal model
      }
      else {
        double m_i, tau_i;
        instream >> m_i >> tau_i;      
        m.push_back(m_i);
        tau.push_back(tau_i);
      }
    }
  }
  int numBck = m.size();

  cout << "n = " << n << endl;
  cout << "s = " << s << endl;
  double btotHat = 0.;
  for (int i=0; i<m.size(); i++){
	btotHat += m[i]/tau[i];
    cout << "m[" << i << "]  = " << m[i] << ",   tau[" << i << "]  = " 
         << tau[i] << endl;
  }
  cout << "estimated total background = " << btotHat << endl;
  cout << endl;

  double ZDisc = getDiscoverySignificance(n, s, m, tau);
  cout << "Discovery significance Z   = " << ZDisc << endl;
  cout << "Corresponding p-value  = " << 1. - TMath::Freq(ZDisc) << endl;
  cout << endl;

  double muTest = 1.0;      // for defining lambda(mu)
  double ZExcl = getExclusionSignificance(muTest, n, s, m, tau);
  cout << "Exclusion significance for mu = " << muTest << " is Z = " << ZExcl << endl;
  cout << "Corresponding p-value = " << 1. - TMath::Freq(ZExcl) << endl;

  return 0;

}
