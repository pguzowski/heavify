/*
 * heavify: create a heavy neutrino flux
 * from a dk2nu beamline simulation
 *
 * author: pawel.guzowski@manchester.ac.uk
 *
 */


#if !defined(__CINT__)
// compiling with root, or externally via gcc

// save this for later
#ifndef ROOT_VERSION
#define __COMPILING_OUTSIDE_ROOT__
#endif

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TGenPhaseSpace.h>
#include <TRandom.h>
#include <TMath.h>

#include <iostream>
#include <vector>
#include <map>

#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"

using namespace std;

const int NYBINS = 500; // number of mass bins (from 0 to 500 MeV)

const int kpdg_nue       =   12;  // extended Geant 53
const int kpdg_nuebar    =  -12;  // extended Geant 52
const int kpdg_numu      =   14;  // extended Geant 56
const int kpdg_numubar   =  -14;  // extended Geant 55

const int kpdg_HNe =      12999; // my own definitions
const int kpdg_HNebar =  -12999; //   (only a single HN/anti_HN
const int kpdg_HNmu =     14999; //   but coupled differently
const int kpdg_HNmubar = -14999; //   to nue/numu)

const int kpdg_eplus     =    -11;  // Geant  5
const int kpdg_eminus    =     11;  // Geant  6
const int kpdg_muplus     =   -13;  // Geant  5
const int kpdg_muminus    =    13;  // Geant  6
const int kpdg_pionplus   =   211;  // Geant  8
const int kpdg_pionminus  =  -211;  // Geant  9
const int kpdg_pi0        =   111;  // Geant  9
const int kpdg_k0long     =   130;  // Geant 10  ( K0=311, K0S=310 )
const int kpdg_k0short    =   310;  // Geant 16
const int kpdg_k0mix      =   311;  
const int kpdg_kaonplus   =   321;  // Geant 11
const int kpdg_kaonminus  =  -321;  // Geant 12
const int kpdg_omegaminus =  3334;  // Geant 24
const int kpdg_omegaplus  = -3334;  // Geant 32

const double kRDET = 100.0;   // set to flux per 100 cm radius


// decay types, from dk2nu ndecay specification
/*
  1  KL0 → νe  + π− + e+
  2  KL0 → ν ̄e + π+ + e−
  3  KL0 → νμ  + π− + μ+
  4  KL0 → ν ̄μ + π+ + μ−
  5  K+  → νμ  + μ+
  6  K+  → νe  + π0 + e+
  7  K+  → νμ  + π0 + μ+
  8  K−  → ν ̄μ + μ−
  9  K−  → ν ̄e + π0 + e−
  10 K−  → ν ̄μ + π0 + μ−
  11 μ+  → ν ̄μ + νe + e+
  12 μ−  → ν  + ν ̄e + e−
  13 π+  → νμ  + μ+
  14 π−  → ν ̄μ + μ−
  999 Other
*/
const map<int, int> mother_pdgs = { 
  { 1,  kpdg_k0long    },
  { 2,  kpdg_k0long    },
  { 3,  kpdg_k0long    },
  { 4,  kpdg_k0long    },
  { 5,  kpdg_kaonplus  },
  { 6,  kpdg_kaonplus  },
  { 7,  kpdg_kaonplus  },
  { 8,  kpdg_kaonminus },
  { 9,  kpdg_kaonminus },
  { 10, kpdg_kaonminus },
  { 11, kpdg_muplus    },
  { 12, kpdg_muminus   },
  { 13, kpdg_pionplus  },
  { 14, kpdg_pionminus }
};
const map<int, vector<int>> daughter_pdgs = {
  { 1,  { kpdg_pionminus, kpdg_eplus    } },
  { 2,  { kpdg_pionplus,  kpdg_eminus   } },
  { 3,  { kpdg_pionminus, kpdg_muplus   } },
  { 4,  { kpdg_pionplus,  kpdg_muminus  } },
  { 5,  { kpdg_muplus                   } },
  { 6,  { kpdg_pi0,       kpdg_eplus    } },
  { 7,  { kpdg_pi0,       kpdg_muplus   } },
  { 8,  { kpdg_muminus                  } },
  { 9,  { kpdg_pi0,       kpdg_eminus   } },
  { 10, { kpdg_pi0,       kpdg_muminus  } },
  { 11, { kpdg_nue,       kpdg_eplus    } },
  { 12, { kpdg_nuebar,    kpdg_eminus   } },
  { 13, { kpdg_muplus                   } },
  { 14, { kpdg_muminus                  } }
};

// Get Q-value for a specific decay type, without heavy neutrino
map<int, double> get_Q_vals(const map<int, vector<int>>& daughters) {
  map<int, double> Q_vals;
  for(auto decay_def : daughters) {
    double Q = TDatabasePDG::Instance()->GetParticle(mother_pdgs.at(decay_def.first))->Mass();
    for(auto dpdg : decay_def.second) {
      Q -=  TDatabasePDG::Instance()->GetParticle(dpdg)->Mass();
    }
    Q_vals[decay_def.first] = Q;
  }
  return Q_vals;
}

// Code based on dk2nu/calcLocationWeights.cxx
map<int, pair<double, double>> calcLocWeights(int decay_type, const map<int, TVector3>& locs, const TLorentzVector& parentP4, const TVector3& decay_pos, const TLorentzVector& heavy_nu_mom, double weight) {
  map<int, pair<double, double>> locWeights;
  for(auto i : locs) {
    const TVector3& pos = i.second - decay_pos;
    double solid_angle = (1. - TMath::Cos(TMath::ATan( kRDET / pos.Mag())))/2.;
    if(parentP4.Vect().Mag() > 0.){
      double costh = parentP4.Vect().Dot(pos) / parentP4.Vect().Mag() / pos.Mag();
      if(costh > 1.) costh = 1.;
      if(costh < -1.) costh = -1.;
      double ang_factor = 1. / ( parentP4.Gamma() * (1. - parentP4.Beta() * costh));
      
      TLorentzVector pp;
      pp.SetVectM(pos, 0.);
      pp.Boost(-parentP4.BoostVector());
      TVector3 cm_dir = (1./ pp.Vect().Mag()) * pp.Vect();
      TLorentzVector nhnm;
      nhnm.SetVectM(heavy_nu_mom.Vect().Mag() * cm_dir, heavy_nu_mom.M());
      nhnm.Boost(parentP4.BoostVector());

      if(decay_type == 11 || decay_type == 12) {
        // TODO: muon decay anisotropy
      }
      locWeights[i.first] = make_pair(nhnm.Vect().Mag(), weight * solid_angle *  ang_factor * ang_factor);
    }
    else {
      locWeights[i.first] = make_pair(heavy_nu_mom.Vect().Mag(), weight * solid_angle);
    }
  }
  return locWeights;
}

// based on Geant4: G4KL3DecayChannel.cc
namespace DalitzKaonParameters {
  // TODO: are these coefficients still valid for heavy neutrinos?
  const map<int, double> pLambda = {
    { 6, 0.0286 }, { 9,  0.0286 }, // k+- -> e+-
    { 7, 0.033  }, { 10, 0.033  }, // k+- -> mu+-
    { 1, 0.0300 }, { 2,  0.0300 }, // K0L -> e+-
    { 3, 0.034  }, { 4,  0.034  }  // K0L -> mu+-
  };
  const map<int, double> pXi0 = {
    { 6, -0.35 }, { 9,  -0.35 }, // k+- -> e+-
    { 7, -0.35 }, { 10, -0.35 }, // k+- -> mu+-
    { 1, -0.11 }, { 2,  -0.11 }, // K0L -> e+-
    { 3, -0.11 }, { 4,  -0.11 }  // K0L -> mu+-
  };
};
double DalitzDensity(const TLorentzVector& p4_k, const TLorentzVector& p4_pi, const TLorentzVector& p4_lep, const TLorentzVector& p4_Hnu, double decay_type) {

  double E_pi = p4_pi.E();
  double E_lep = p4_lep.E();
  double E_Hnu = p4_Hnu.E();

  double massK = p4_k.M();
  double massPi = p4_pi.M();
  double massLep = p4_lep.M();
  double massHnu = p4_Hnu.M();

  double E_pi_max = 0.5*(massK * massK + massPi * massPi - massLep * massLep - massHnu * massHnu)/massK;
  double E = E_pi_max - E_pi;
  double Q2 = massK * massK + massPi * massPi - 2. * massK * E_pi;

  const double pLambda = DalitzKaonParameters::pLambda.at(decay_type);
  const double pXi0 = DalitzKaonParameters::pXi0.at(decay_type);


  double F = 1. + pLambda * Q2 / massPi / massPi;
  //double Fmax = 1.;
  //if(pLambda > 0.)
  double Fmax = (1. + pLambda * (massK * massK / massPi / massPi + 1.));

  double Xi = pXi0 * (1. + pLambda * Q2 / massPi / massPi);

  double coeffA = massK * (2. * E_lep * E_Hnu - massK * E) + massLep * massLep * (E/4. - E_Hnu);
  double coeffB = massLep * massLep * (E_Hnu - E / 2.);
  double coeffC = massLep * massLep * E / 4.;

  double RhoMax = (Fmax * Fmax) * (massK * massK * massK / 8.);
  double Rho = (F * F) * (coeffA + coeffB * Xi + coeffC * Xi * Xi);

  return Rho / RhoMax;
}

TLorentzVector make_decay(TGenPhaseSpace* decayer, const TLorentzVector& parentP4, const map<int, TVector3>& locations, map<int, pair<double, double>>& locWeights, const TVector3& decay_pos, int decay_type) {

  // Generate uniformly [This can infinite loop]
  const double wmax = decayer->GetWtMax();
  double weight = -1.;
  const unsigned int nmaxthrows = 100000;
  unsigned int nthrows = 0;
  while(gRandom->Uniform(wmax) > weight) {
    weight = decayer->Generate();

    // Kaon dalitz decay
    double ddw_weight = weight;
    if(decay_type <= 4 || decay_type == 6 || decay_type == 7 || decay_type == 9 || decay_type == 10 ) {
      double ddw = DalitzDensity(parentP4, *(decayer->GetDecay(1)), *(decayer->GetDecay(2)), *(decayer->GetDecay(0)), decay_type);
      if(gRandom->Uniform() > ddw) {
        weight = -1.; // will cause rethrow at while(gRandom->Uniform(wmax) > weight)
        ddw_weight *= ddw;
      }
    }

    // muon decay
    if(decay_type == 11 || decay_type == 12) {
      // TODO: Muon decay beta-style energy spectrum
    }
    
    // break infintie loop
    if(nthrows++ > nmaxthrows) {
      if(weight < 0.) {
        // kaon dalitz decay
        weight = ddw_weight;
      }
      break;
    }
  }

  TLorentzVector nuP4 = *(decayer->GetDecay(0));
  if(nthrows >= nmaxthrows) {
    locWeights = calcLocWeights(decay_type, locations, parentP4, decay_pos, nuP4, weight/wmax);
    nuP4.Boost(parentP4.BoostVector());
    locWeights[0] = make_pair(nuP4.Vect().Mag(), weight / wmax);
  }
  else {
    locWeights = calcLocWeights(decay_type, locations, parentP4, decay_pos, nuP4, 1.);
    nuP4.Boost(parentP4.BoostVector());
    locWeights[0] = make_pair(nuP4.Vect().Mag(), 1.);
  }
  return nuP4;
}


map<int, map<double, TGenPhaseSpace*>> decayers;

// create decayer if it doesn't exist
TGenPhaseSpace* get_decayer(int decay_type, double heavy_mass) {
  if(daughter_pdgs.find(decay_type) == daughter_pdgs.end()) {
    throw "Decayer failed";
  }
  map<double, TGenPhaseSpace*>& decay = decayers[decay_type];
  if(decay.find(heavy_mass) != decay.end()) {
    return decay[heavy_mass];
  }

  vector<double> daughter_masses = { heavy_mass };
  for(auto dpdg : daughter_pdgs.at(decay_type)) {
    daughter_masses.push_back(TDatabasePDG::Instance()->GetParticle(dpdg)->Mass());
  }

  TGenPhaseSpace *decayer = new TGenPhaseSpace;
  TLorentzVector dummy(0.,0.,0.,TDatabasePDG::Instance()->GetParticle(mother_pdgs.at(decay_type))->Mass());
  decayer->SetDecay(dummy, daughter_masses.size(), daughter_masses.data());
  decay[heavy_mass] = decayer;
  return decayer;
}

// prepare histograms if they don't exist
map<int, TH1*>& GetHist(map<int, map<int, map<int, TH1*>>>& hists, int nupdg, int ndecay, const map<int, string>& loc_names) {
  
  if(hists[nupdg].find(ndecay) == hists[nupdg].end()) {

    const map<int, string> nunames = { { kpdg_numubar, "anumu" }, { kpdg_nuebar, "anue" }, { kpdg_nue, "nue" }, { kpdg_numu, "numu" },
                                       { kpdg_HNmubar, "aHNmu" }, { kpdg_HNebar, "aHNe" }, { kpdg_HNe, "HNe" }, { kpdg_HNmu, "HNmu" } };
  
    if(abs(nupdg) > 100) {
      hists[nupdg][ndecay][0] = new TH2D(Form("h_%s_d%02d",nunames.at(nupdg).c_str(),ndecay), ";Heavy Neutrino Momentum (GeV/c);Heavy Neutrino Mass (GeV/c^{2})", 1000, 0, 10, NYBINS, 0., 0.5);
    }
    else {
      hists[nupdg][ndecay][0] = new TH1D(Form("h_%s_d%02d",nunames.at(nupdg).c_str(),ndecay), ";Light Neutrino Momentum (GeV/c);", 1000, 0, 10);
    }
    for(auto l : loc_names) {
      if(abs(nupdg) > 100) {
        hists[nupdg][ndecay][l.first] = new TH2D(Form("h_%s_d%02d_%s", nunames.at(nupdg).c_str(),ndecay, l.second.c_str()), Form("%s;Heavy Neutrino Momentum (GeV/c);Heavy Neutrino Mass (GeV/c^{2})",l.second.c_str()), 100, 0, 10, NYBINS, 0, 0.5);
      }
      else {
        hists[nupdg][ndecay][l.first] = new TH1D(Form("h_%s_d%02d_%s", nunames.at(nupdg).c_str(),ndecay, l.second.c_str()), Form("%s;Light Neutrino Momentum (GeV/c)",l.second.c_str()), 100, 0, 10);
      }
    }
    for(auto h : hists[nupdg][ndecay]) h.second->Sumw2();


  }
  return hists[nupdg][ndecay];
}


/*****
 * MAIN FUNCTION
 *****/

void heavify(const char* infilename = 0, const char* outfilename = 0) {

  if (!infilename || !outfilename) {
    cout << "Usage: heavify(input_filename, output_filename);" << endl;
    // emulate root, so as not to confuse users
    cout << "root [1] " << flush;
    return;
  }

  TFile *f = new TFile(infilename);
  if(!f || f->IsZombie()) {
    cerr << "Filename " << infilename << " does not exist!" << endl;
    return;
  }


  TTree *mt = (TTree*)f->Get("dkmetaTree");
  if(!mt) {
    cerr << "File " << infilename << " does not contain a dkmetaTree!" << endl;
    return;
  }
  bsim::DkMeta *met = 0;
  mt->SetBranchAddress("dkmeta", &met);
  mt->GetEntry(0);

  map<int, TVector3> locations;
  map<int, string> loc_names;
  {
    int FL = 1;
    for(auto l : met->location) {
      locations[FL] = TVector3(l.x,l.y,l.z);
      loc_names[FL] = l.name;
      ++FL;
    }
  }


  TTree *t = (TTree*)f->Get("dk2nuTree");
  if(!t) {
    cerr << "File " << infilename << " does not contain a dk2nuTree!" << endl;
    return;
  }
  bsim::Dk2Nu *dk = 0;
  t->SetBranchAddress("dk2nu", &dk);

  const map<int, vector<int>>& dpdgs = daughter_pdgs;
  const map<int, double> Qv = get_Q_vals(dpdgs);

  map<int, map<int, map<int, TH1*>>> hists;
  map<int, TH1*> hists_ndecay;

  TFile *fout = new TFile(outfilename,"recreate");

  map<int, string> nus = { {kpdg_numubar, "anumu"}, {kpdg_nuebar, "anue"}, {kpdg_nue, "nue"}, {kpdg_numu, "numu"} };
  for(auto p : nus) {
    hists_ndecay[p.first] = new TH2D(Form("h_ndecay_%s",p.second.c_str()), "", 16, 0, 16, NYBINS, 0, 0.500);
    hists_ndecay[p.first]->Sumw2();
  }



  for (int i = 0; i < t->GetEntries(); ++i){
    if((i+1)%10000 == 0) {
      cout << "Processed " << (i+1) << " entries" << endl;
    }
    t->GetEntry(i);
    const int nu_pdg = dk->decay.ntype;
    const int decay_type = dk->decay.ndecay;
    const int parent_pdg = dk->decay.ptype;
    const double parent_mass = TDatabasePDG::Instance()->GetParticle(parent_pdg)->Mass();
    if(Qv.find(decay_type) == Qv.end()) {
      cout << "Unknown decay: " << decay_type << endl;
      continue;
    }
    TVector3 decay_pos(dk->decay.vx, dk->decay.vy, dk->decay.vz);

    const int HN_pdg = 1000*nu_pdg + (nu_pdg>0?999:-999);

    map<int, TH1*>& h_nu = GetHist(hists, nu_pdg, decay_type, loc_names);
    map<int, TH1*>& h_HN = GetHist(hists, HN_pdg, decay_type, loc_names);
    
    for(unsigned int k = 0; k <= locations.size(); ++k) {
      h_nu[k]->Fill(dk->nuray[k].E, dk->nuray[k].wgt);
    }
    for(int M = 1; M <= NYBINS; ++M) {
      const double heavy_mass = h_HN[0]->GetYaxis()->GetBinCenter(M);
      if(Qv.at(decay_type) > heavy_mass) {
        hists_ndecay[nu_pdg]->Fill(decay_type, heavy_mass);

        TLorentzVector parentP4;
        parentP4.SetXYZM(dk->decay.pdpx, dk->decay.pdpy, dk->decay.pdpz, parent_mass);
        map<int, pair<double, double>> locWeights;
        TLorentzVector hnu = make_decay(get_decayer(decay_type, heavy_mass), parentP4, locations, locWeights, decay_pos, decay_type);

        for(auto l : locWeights) {
          dynamic_cast<TH2*>(h_HN[l.first])->Fill(l.second.first, heavy_mass, l.second.second);
        }
      }
    }
  }
  for(auto i : hists) for(auto j : i.second) for(auto k : j.second) {
    k.second->Write();
  }
  for(auto i : hists_ndecay) i.second->Write();
  fout->Close();
}

// finally add a main() function if compiling outside root
// usage: ./heavify -i input_file -o output_file
#ifdef __COMPILING_OUTSIDE_ROOT__
#include <unistd.h>
int main(int argc, char** argv) {
  string infname, outfname;
  char c;
  while((c = getopt(argc, argv, "i:o:h" )) != -1) {
    switch(c) {
      case 'i':
        infname = optarg;
        break;
      case 'o':
        outfname = optarg;
        break;
      case 'h':
        cout << "Usage: " << endl;
        cout << "\t ./heavify -i <input file> -o <output file>" << endl;
        return -1;
      default:
        break;
    }
  }
  if(!infname.empty() && !outfname.empty()) {
    heavify(infname.c_str(), outfname.c_str());
  }
  return 0;
}
#endif

#else // if __CINT__ is defined
// Running inside root interpreter

void heavify(const char* = 0, const char* = 0) {

  // Assumes dk2nu has been loaded via ups
  if(gSystem->Getenv("DK2NU") == 0) {
    cerr << "Need to setup dk2nu via ups first!" << endl;
    return;
  }
  string dk2nu_lib = gSystem->Getenv("DK2NU_LIB");
  string incdir = " -I";
  string dk2nu_inc = gSystem->Getenv("DK2NU_INC");
  gSystem->AddIncludePath((incdir+dk2nu_inc+" ").c_str());
  gSystem->Load((dk2nu_lib+"/libdk2nuTree.so").c_str());
  gSystem->Load("libEG");

  // reload this file, with compilation
  // need to have it on a timer, or root complains that
  // it cannot unload this same file before
  // recompiling
  (new TTimer("gROOT->LoadMacro(\"heavify.C+\"); heavify();", -1))->Start(-1,kTRUE);
}

#endif
