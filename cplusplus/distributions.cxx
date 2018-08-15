#include <ROOT/RDataFrame.hxx>

#include <TInterpreter.h>

#include "common_definitions.h"
#include "parameters_global.h"
#include "common_algorithms.h"
#include "parameters.h"
#include "common.h"

#include <iostream>

#include <TROOT.h>

using RDF = ROOT::RDataFrame;
using namespace std;

// Get input (temporal options)
auto treename          = "distilled";
auto selected_diagonal = "d45b_56t";
auto prefix            = "";
auto outputDir         = ".";

int main(int argc, char **argv)
{
  gInterpreter->Declare(R"cpp(
    #include "common_definitions.h"
    #include "parameters_global.h"
    #include "common_algorithms.h"
    #include "parameters.h"
    #include "common.h"
  )cpp");

  // Read input files
  auto fname = argv[1];
  auto input_file  = fname ; // Created with distill.py
  gInterpreter->ProcessLine("Init(\"45b_56t\")");

  // default parameters
  unsigned int detailsLevel = 0;     // 0: no details, 1: some details, >= 2 all details
  bool overrideCutSelection = false;    // whether the default cut selection should be overriden by the command-line selection
  string cutSelectionString;
  string outputDir = ".";
  string inputDir = ".";
  double input_n_si = 4.0;
  int time_group_divisor = 0;
  int time_group_remainder = 0;
  int event_group_divisor = 0;
  int event_group_index = 0;
  unsigned int evIdxStep = 1;
  unsigned int maxTaggedEvents = 0;    // 0 means no maximum

  printf("* detailsLevel = %u\n", detailsLevel);
  printf("* outputDir = %s\n", outputDir.c_str());
  printf("* inputDir = %s\n", inputDir.c_str());
  printf("* input n_si = %.3f\n", input_n_si);
  printf("* time_group_divisor = %i\n", time_group_divisor);
  printf("* time_group_remainder = %i\n", time_group_remainder);
  printf("* event_group_divisor = %i\n", event_group_divisor);
  printf("* event_group_index = %i\n", event_group_index);
  printf("* evIdxStep = %u\n", evIdxStep);
  printf("* maxTaggedEvents = %u\n", maxTaggedEvents);

  // select cuts
  auto anal = (Analysis*) gInterpreter->ProcessLine("anal;");
  anal->BuildCuts();
  anal->n_si = input_n_si;

  // print info
  printf("\n");
  printf("------------------------------ environment ------------------------------\n");
  env.Print();
  printf("\n");
  printf("------------------------------- analysis --------------------------------\n");
  anal->Print();
  printf("\n");

  // alignment init
  for (unsigned int i = 0; i < alignmentSources.size(); ++i)
  {
    printf("\n---------- alignment source %u ----------\n", i);
    alignmentSources[i].Init();
  }
  printf("\n\n");

  // binnings
  vector<string> binnings;
  binnings.push_back("ub");
  //binnings.push_back("eb");
  binnings.push_back("ob-1-10-0.2");
  binnings.push_back("ob-1-30-0.2");


  //#########################################################
  //######   READ INPUT FILE, INITIALIZE RDATAFRAME    ######
  //#########################################################

  // Read all branches
  RDF rdf(treename, input_file);

  // book metadata histograms
  unsigned int timestamp_bins = timestamp_max - timestamp_min + 1.;

  Binning binning_setup[binnings.size()];
  map<unsigned int, ROOT::RDF::RResultPtr<TH1D>> bh_t_Nev_before, bh_t_Nev_after_no_corr;
	map<unsigned int, ROOT::RDF::RResultPtr<TH1D>> bh_t_before, bh_t_after_no_corr, bh_t_after;
	//map<unsigned int, TProfile*> bp_t_phi_corr, bp_t_full_corr;

  for (unsigned int bi = 0; bi < binnings.size(); ++bi)
  {
    Binning b;
    BuildBinningRDF(*anal, binnings[bi], b);
    binning_setup[bi] = b;
  }
  // zero counters
     unsigned long n_ev_full = 0;
     map<unsigned int, unsigned long> n_ev_cut;
     for (unsigned int ci = 1; ci <= anal->N_cuts; ++ci)
         n_ev_cut[ci] = 0;

     double th_min = 1E100;
     double th_y_L_min = +1E100, th_y_R_min = +1E100;

     unsigned int N_anal=0, N_anal_zeroBias=0;
     unsigned int N_zeroBias_el=0, N_zeroBias_el_RP_trig=0;
     unsigned int N_4outof4=0, N_el=0;
     unsigned int N_el_T2trig=0, N_4outof4_T2trig=0;
     unsigned int N_el_raw=0;

  // #########################################################
  // ###### FILTER, BUILD HISTOGRAMS - START EVENT LOOP ######
  // #########################################################


  auto f1 = rdf.Filter("! SkipTime( timestamp )", "check time - selected");

  // Diagonal cut (L831)
  auto f2 = f1.Filter("v_L_2_F && v_L_2_N && v_R_2_F && v_R_2_N", "allDiagonalRPs");

  auto model = ROOT::RDF::TH1DModel("h_timestamp_dgn", ";timestamp;rate   (Hz)", int(timestamp_bins), timestamp_min-0.5, timestamp_max+0.5);
  auto h_timestamp_dgn = f2.Histo1D(model, "timestamp");

  // Not cut for this filter in original code
  auto f_zerobias = f2.Filter("! ((trigger_bits & 512) != 0)", "zero_bias_event");

  auto r2 = f2.Define("h_al", "ApplyFineAlignment( timestamp , x_L_1_F, x_L_2_N, x_L_2_F, x_R_1_F, x_R_2_N, x_R_2_F, y_L_1_F, y_L_2_N, y_L_2_F, y_R_1_F, y_R_2_N, y_R_2_F)" )
               .Define("h_al_x_L_1_F", "h_al.L_1_F.x")
               .Define("h_al_x_L_2_N", "h_al.L_2_N.x")
               .Define("h_al_x_L_2_F", "h_al.L_2_F.x")
               .Define("h_al_y_L_1_F", "h_al.L_1_F.y")
               .Define("h_al_y_L_2_N", "h_al.L_2_N.y")
               .Define("h_al_y_L_2_F", "h_al.L_2_F.y")
               .Define("h_al_x_R_1_F", "h_al.R_1_F.x")
               .Define("h_al_x_R_2_N", "h_al.R_2_N.x")
               .Define("h_al_x_R_2_F", "h_al.R_2_F.x")
               .Define("h_al_y_R_1_F", "h_al.R_1_F.y")
               .Define("h_al_y_R_2_N", "h_al.R_2_N.y")
               .Define("h_al_y_R_2_F", "h_al.R_2_F.y");

   // fill pre-selection histograms (Line 860 - 866)
   map<int, ROOT::RDF::TH2DModel>al_nosel_models = {
      {0, ROOT::RDF::TH2DModel("h_y_L_1_F_vs_x_L_1_F_al_nosel", ";x^{L,1,F};y^{L,1,F}", 150, -15., 15., 300, -30., +30.)},
      {1 , ROOT::RDF::TH2DModel("h_y_L_2_N_vs_x_L_2_N_al_nosel", ";x^{L,2,N};y^{L,2,N}", 150, -15., 15., 300, -30., +30.)},
      {2 , ROOT::RDF::TH2DModel("h_y_L_2_F_vs_x_L_2_F_al_nosel", ";x^{L,2,F};y^{L,2,F}", 150, -15., 15., 300, -30., +30.)},
      {3 , ROOT::RDF::TH2DModel("h_y_R_1_F_vs_x_R_1_F_al_nosel", ";x^{R,1,F};y^{R,1,F}", 150, -15., 15., 300, -30., +30.)},
      {4 , ROOT::RDF::TH2DModel("h_y_R_2_N_vs_x_R_2_N_al_nosel", ";x^{R,2,N};y^{R,2,N}", 150, -15., 15., 300, -30., +30.)},
      {5 , ROOT::RDF::TH2DModel("h_y_R_2_F_vs_x_R_2_F_al_nosel", ";x^{R,2,F};y^{R,2,F}", 150, -15., 15., 300, -30., +30.)}
   };

   auto h_y_L_1_F_vs_x_L_1_F_al_nosel = r2.Histo2D(al_nosel_models[0], "h_al_x_L_1_F", "h_al_y_L_1_F");
   auto h_y_L_2_N_vs_x_L_2_N_al_nosel = r2.Histo2D(al_nosel_models[1], "h_al_x_L_2_N", "h_al_y_L_2_N");
   auto h_y_L_2_F_vs_x_L_2_F_al_nosel = r2.Histo2D(al_nosel_models[2], "h_al_x_L_2_F", "h_al_y_L_2_F");
   auto h_y_R_1_F_vs_x_R_1_F_al_nosel = r2.Histo2D(al_nosel_models[3], "h_al_x_R_1_F", "h_al_y_R_1_F");
   auto h_y_R_2_N_vs_x_R_2_N_al_nosel = r2.Histo2D(al_nosel_models[4], "h_al_x_R_2_N", "h_al_y_R_2_N");
   auto h_y_R_2_F_vs_x_R_2_F_al_nosel = r2.Histo2D(al_nosel_models[5], "h_al_x_R_2_F", "h_al_y_R_2_F");

   auto r3 = r2.Define("kinematics", "DoReconstruction( h_al )");

   auto r4 = r3.Define("k_th_x_R",        "kinematics.th_x_R")
          .Define("k_th_y_R",        "kinematics.th_y_R")
          .Define("k_th_x_L",        "kinematics.th_x_L")
          .Define("k_th_y_L",        "kinematics.th_y_L")
          .Define("k_th_x",          "kinematics.th_x")
          .Define("k_th_y",          "kinematics.th_y")
          .Define("minus_k_th_y",    "- kinematics.th_y")
          .Define("k_vtx_x",         "kinematics.vtx_x")
          .Define("k_vtx_x_L",       "kinematics.vtx_x_L")
          .Define("k_vtx_x_R",       "kinematics.vtx_x_R")
          .Define("k_vtx_y",         "kinematics.vtx_y")
          .Define("k_vtx_y_L",       "kinematics.vtx_y_L")
          .Define("k_vtx_y_R",       "kinematics.vtx_y_R")
          .Define("k_th_y_L_F",      "kinematics.th_y_L_F")
          .Define("k_th_y_L_N",      "kinematics.th_y_L_N")
          .Define("k_th_y_R_F",      "kinematics.th_y_R_F")
          .Define("k_th_y_R_N",      "kinematics.th_y_R_N")
          .Define("k_th_x_diffLR",   "kinematics.th_x_R - kinematics.th_x_L")
          .Define("k_th_y_diffLR",   "kinematics.th_y_R - kinematics.th_y_L")
          .Define("k_th_x_diffLF",   "kinematics.th_x_L - kinematics.th_x")
          .Define("k_th_x_diffRF",   "kinematics.th_x_R - kinematics.th_x")
          .Define("k_th_y_L_diffNF", "kinematics.th_y_L_F - kinematics.th_y_L_N")
          .Define("k_th_y_R_diffNF", "kinematics.th_y_R_F - kinematics.th_y_R_N")
          .Define("k_vtx_x_diffLR",  "kinematics.vtx_x_R - kinematics.vtx_x_L")
          .Define("k_vtx_y_diffLR",  "kinematics.vtx_y_R - kinematics.vtx_y_L")
          .Define("k_t",             "kinematics.t")
          .Define("k_th",            "kinematics.th")
          .Define("k_phi",           "kinematics.phi");

    auto r5 = r4.Define("cutdata", "EvaluateCutsRDF( h_al, kinematics )");

    // Elastic cut
    auto f4 = r5.Filter("cutdata.select", "elastic cut");

    // Define normalization and norm_corr colums
    auto r6 = r5.Define("norm_corr",     "getNorm_corr( timestamp )" )
           .Define("normalization", "getNormalization( norm_corr )");

    auto h_timestamp_sel = f4.Histo1D(ROOT::RDF::TH1DModel("h_timestamp_sel", ";timestamp;rate   (Hz)", int(timestamp_bins), timestamp_min-0.5, timestamp_max+0.5), "timestamp");

    // fill histograms
    map<int, ROOT::RDF::TH2DModel> noal_sel_models = {
        {0, ROOT::RDF::TH2DModel("h_y_L_1_F_vs_x_L_1_F_noal_sel", ";x^{L,1,F};y^{L,1,F}", 100, -3., +3., 300, -30., +30.)},
        {1, ROOT::RDF::TH2DModel("h_y_L_2_N_vs_x_L_2_N_noal_sel", ";x^{L,2,N};y^{L,2,N}", 100, -3., +3., 300, -30., +30.)},
        {2, ROOT::RDF::TH2DModel("h_y_L_2_F_vs_x_L_2_F_noal_sel", ";x^{L,2,F};y^{L,2,F}", 100, -3., +3., 300, -30., +30.)},
        {3, ROOT::RDF::TH2DModel("h_y_R_1_F_vs_x_R_1_F_noal_sel", ";x^{R,1,F};y^{R,1,F}", 100, -3., +3., 300, -30., +30.)},
        {4, ROOT::RDF::TH2DModel("h_y_R_2_N_vs_x_R_2_N_noal_sel", ";x^{R,2,N};y^{R,2,N}", 100, -3., +3., 300, -30., +30.)},
        {5, ROOT::RDF::TH2DModel("h_y_R_2_F_vs_x_R_2_F_noal_sel", ";x^{R,2,F};y^{R,2,F}", 100, -3., +3., 300, -30., +30.)}
    };

    auto h_y_L_1_F_vs_x_L_1_F_noal_sel = f4.Histo2D(noal_sel_models[0], "x_L_1_F", "y_L_1_F");
    auto h_y_L_2_N_vs_x_L_2_N_noal_sel = f4.Histo2D(noal_sel_models[1], "x_L_2_N", "y_L_2_N");
    auto h_y_L_2_F_vs_x_L_2_F_noal_sel = f4.Histo2D(noal_sel_models[2], "x_L_2_F", "y_L_2_F");
    auto h_y_R_1_F_vs_x_R_1_F_noal_sel = f4.Histo2D(noal_sel_models[3], "x_R_1_F", "y_R_1_F");
    auto h_y_R_2_N_vs_x_R_2_N_noal_sel = f4.Histo2D(noal_sel_models[4], "x_R_2_N", "y_R_2_N");
    auto h_y_R_2_F_vs_x_R_2_F_noal_sel = f4.Histo2D(noal_sel_models[5], "x_R_2_F", "y_R_2_F");

    map<int, ROOT::RDF::TH2DModel> al_sel_models = {
        {0, ROOT::RDF::TH2DModel("h_y_L_1_F_vs_x_L_1_F_al_sel", ";x^{L,1,F};y^{L,1,F}", 100, -3., +3., 300, -30., +30.)},
        {1, ROOT::RDF::TH2DModel("h_y_L_2_N_vs_x_L_2_N_al_sel", ";x^{L,2,N};y^{L,2,N}", 100, -3., +3., 300, -30., +30.)},
        {2, ROOT::RDF::TH2DModel("h_y_L_2_F_vs_x_L_2_F_al_sel", ";x^{L,2,F};y^{L,2,F}", 100, -3., +3., 300, -30., +30.)},
        {3, ROOT::RDF::TH2DModel("h_y_R_1_F_vs_x_R_1_F_al_sel", ";x^{R,1,F};y^{R,1,F}", 100, -3., +3., 300, -30., +30.)},
        {4, ROOT::RDF::TH2DModel("h_y_R_2_N_vs_x_R_2_N_al_sel", ";x^{R,2,N};y^{R,2,N}", 100, -3., +3., 300, -30., +30.)},
        {5, ROOT::RDF::TH2DModel("h_y_R_2_F_vs_x_R_2_F_al_sel", ";x^{R,2,F};y^{R,2,F}", 100, -3., +3., 300, -30., +30.)}
    };

    auto h_y_L_1_F_vs_x_L_1_F_al_sel = f4.Histo2D(al_sel_models[0], "h_al_x_L_1_F", "h_al_y_L_1_F");
    auto h_y_L_2_N_vs_x_L_2_N_al_sel = f4.Histo2D(al_sel_models[1], "h_al_x_L_2_N", "h_al_y_L_2_N");
    auto h_y_L_2_F_vs_x_L_2_F_al_sel = f4.Histo2D(al_sel_models[2], "h_al_x_L_2_F", "h_al_y_L_2_F");
    auto h_y_R_1_F_vs_x_R_1_F_al_sel = f4.Histo2D(al_sel_models[3], "h_al_x_R_1_F", "h_al_y_R_1_F");
    auto h_y_R_2_N_vs_x_R_2_N_al_sel = f4.Histo2D(al_sel_models[4], "h_al_x_R_2_N", "h_al_y_R_2_N");
    auto h_y_R_2_F_vs_x_R_2_F_al_sel = f4.Histo2D(al_sel_models[5], "h_al_x_R_2_F", "h_al_y_R_2_F");

    // Line 1157 (k.th_x_R - k.th_x_L)
    //           (k.th_y_R - k.th_y_L)
    auto th_x_diffLR = f4.Histo1D(ROOT::RDF::TH1DModel("th_x_diffLR", ";#theta_{x}^{R} - #theta_{x}^{L}", 1000, -500E-6, +500E-6), "k_th_x_diffLR");
    auto th_y_diffLR = f4.Histo1D(ROOT::RDF::TH1DModel("th_y_diffLR", ";#theta_{y}^{R} - #theta_{y}^{L}", 500, -50E-6, +50E-6), "k_th_y_diffLR");

    // Line 1160 (k.th_x_L - k.th_x)
    //           (k.th_x_R - k.th_x)
    auto th_x_diffLF = f4.Histo1D(ROOT::RDF::TH1DModel("th_x_diffLF", ";#theta_{x}^{L} - #theta_{x}", 400, -200E-6, +200E-6), "k_th_x_diffLF");
    auto th_x_diffRF = f4.Histo1D(ROOT::RDF::TH1DModel("th_x_diffRF", ";#theta_{x}^{R} - #theta_{x}", 400, -200E-6, +200E-6), "k_th_x_diffRF");

    // Line 1163 (k.th_x, k.th_x_R - k.th_x_L)
    //           (k.th_y, k.th_y_R - k.th_y_L)
    //           (k.vtx_x, k.th_x_R - k.th_x_L)
    auto h_th_x_diffLR_vs_th_x  = f4.Histo2D(ROOT::RDF::TH2DModel("h_th_x_diffLR_vs_th_x", ";#theta_{x};#theta_{x}^{R} - #theta_{x}^{L}", 100, -300E-6, +300E-6, 120, -120E-6, +120E-6), "k_th_x", "k_th_x_diffLR");
    auto h_th_y_diffLR_vs_th_y  = f4.Histo2D(ROOT::RDF::TH2DModel("h_th_y_diffLR_vs_th_y", ";#theta_{y};#theta_{y}^{R} - #theta_{y}^{L}", 100, -500E-6, +500E-6, 120, -120E-6, +120E-6), "k_th_y", "k_th_y_diffLR");
    auto h_th_x_diffLR_vs_vtx_x = f4.Histo2D(ROOT::RDF::TH2DModel("h_th_x_diffLR_vs_vtx_x", ";vtx_{x};#theta_{x}^{R} - #theta_{x}^{L}", 100, -300E-3, +300E-3, 120, -120E-6, +120E-6), "k_vtx_x", "k_th_x_diffLR");

    auto h_th_y_L_vs_th_x_L = f4.Histo2D(ROOT::RDF::TH2DModel("h_th_y_L_vs_th_x_L", ";#theta_{x}^{L};#theta_{y}^{L}", 100, -115E-6, +11E-5, 100, 22E-6, +102E-6), "k_th_x_L", "k_th_y_L");
    auto h_th_y_R_vs_th_x_R = f4.Histo2D(ROOT::RDF::TH2DModel("h_th_y_R_vs_th_x_R", ";#theta_{x}^{R};#theta_{y}^{R}", 100, -125E-6, +12E-5, 100, 27E-6, +102E-6), "k_th_x_R", "k_th_y_R");
    auto h_th_y_vs_th_x     = f4.Histo2D(ROOT::RDF::TH2DModel("h_th_y_vs_th_x", ";#theta_{x};#theta_{y}", 100, -300E-6, +300E-6, 100, -150E-6, +150E-6), "k_th_x", "k_th_y");

    auto h_th_y_L_vs_th_y_R = f4.Histo2D(ROOT::RDF::TH2DModel("h_th_y_L_vs_th_y_R", ";#theta_{y}^{R};#theta_{y}^{L}",
                                     300, -150E-6, +150E-6, 300, -150E-6, +150E-6), "k_th_y_R", "k_th_y_L");

    // Line 1203: (k.th_x)
    //            (k.th_y)
    auto h_th_x = f4.Histo1D(ROOT::RDF::TH1DModel("h_th_x", ";#theta_{x}", 250, -500E-6, +500E-6), "k_th_x");
    auto h_th_y = f4.Histo1D(ROOT::RDF::TH1DModel("h_th_y", ";#theta_{y}", 250, -500E-6, +500E-6), "k_th_y");

    // Line 1205: (-k.th_y)
    auto h_th_y_flipped = f4.Histo1D(ROOT::RDF::TH1DModel("h_th_y_flipped", ";#theta_{y}", 250, -500E-6, +500E-6), "minus_k_th_y");

    // Line 1207: (k.th_x_L)
    //            (k.th_x_R)
    auto h_th_x_L = f4.Histo1D(ROOT::RDF::TH1DModel("h_th_x_L", ";#theta_{x}^{L}", 250, -500E-6, +500E-6), "k_th_x_L");
    auto h_th_x_R = f4.Histo1D(ROOT::RDF::TH1DModel("h_th_x_R", ";#theta_{x}^{R}", 250, -500E-6, +500E-6), "k_th_x_R");

    auto h_th_y_L = f4.Histo1D(ROOT::RDF::TH1DModel("h_th_y_L", ";#theta_{y}^{L}", 250, -500E-6, +500E-6), "k_th_y_L");
    auto h_th_y_R = f4.Histo1D(ROOT::RDF::TH1DModel("h_th_y_R", ";#theta_{y}^{R}", 250, -500E-6, +500E-6), "k_th_y_R");

    // Line 1213: (k.th_y_L_F)
    //            (k.th_y_L_N)
    //            (k.th_y_R_N)
    //            (k.th_y_R_F)
    auto h_th_y_L_F = f4.Histo1D(ROOT::RDF::TH1DModel("h_th_y_L_F", ";#theta_{y}^{L_F}", 250, -500E-6, +500E-6), "k_th_y_L_F");
    auto h_th_y_L_N = f4.Histo1D(ROOT::RDF::TH1DModel("h_th_y_L_N", ";#theta_{y}^{L_N}", 250, -500E-6, +500E-6), "k_th_y_L_N");
    auto h_th_y_R_N = f4.Histo1D(ROOT::RDF::TH1DModel("h_th_y_R_N", ";#theta_{y}^{R_N}", 250, -500E-6, +500E-6), "k_th_y_R_N");
    auto h_th_y_R_F = f4.Histo1D(ROOT::RDF::TH1DModel("h_th_y_R_F", ";#theta_{y}^{R_F}", 250, -500E-6, +500E-6), "k_th_y_R_F");

    // fill vertex histograms

    // Line 1220 (k.vtx_x)
    //           (k.vtx_x_L)
    //           (k.vtx_x_R)
    auto h_vtx_x    = f4.Histo1D(ROOT::RDF::TH1DModel("h_vtx_x", ";x^{*}", 100, -0.5, +0.5) , "k_vtx_x");
    auto h_vtx_x_L  = f4.Histo1D(ROOT::RDF::TH1DModel("h_vtx_x_L", ";x^{*,L}", 100, -0.5, +0.5) , "k_vtx_x_L");
    auto h_vtx_x_R  = f4.Histo1D(ROOT::RDF::TH1DModel("h_vtx_x_R", ";x^{*,R}", 100, -0.5, +0.5) , "k_vtx_x_R");

    // Line 1224 (k.vtx_y)
    //           (k.vtx_y_L)
    //           (k.vtx_y_R)
    auto h_vtx_y    = f4.Histo1D(ROOT::RDF::TH1DModel("h_vtx_y", ";y^{*}", 100, -0.5, +0.5), "k_vtx_y");
    auto h_vtx_y_L  = f4.Histo1D(ROOT::RDF::TH1DModel("h_vtx_y_L", ";y^{*,L}", 100, -0.5, +0.5), "k_vtx_y_L");
    auto h_vtx_y_R  = f4.Histo1D(ROOT::RDF::TH1DModel("h_vtx_y_R", ";y^{*,R}", 100, -0.5, +0.5), "k_vtx_y_R");

    // Line 1228:
    //            (k.vtx_x_R, k.vtx_x_L)
    //            (k.vtx_y_R, k.vtx_y_L)
    auto h_vtx_x_L_vs_vtx_x_R = f4.Histo2D(ROOT::RDF::TH2DModel("h_vtx_x_L_vs_vtx_x_R", ";x^{*,R};x^{*,L}", 100, -0.5, +0.5, 100, -0.5, +0.5), "k_vtx_x_R", "k_vtx_x_L");
    auto h_vtx_y_L_vs_vtx_y_R = f4.Histo2D(ROOT::RDF::TH2DModel("h_vtx_y_L_vs_vtx_y_R", ";y^{*,R};y^{*,L}", 100, -0.5, +0.5, 100, -0.5, +0.5), "k_vtx_y_R", "k_vtx_y_L");

    // Line 1231:
    //            (k.th_x_L, k.vtx_x_L)
    //            (k.th_x_R, k.vtx_x_R)
    //            (k.th_y_L, k.vtx_y_L)
    //            (k.th_y_R, k.vtx_y_R)
    auto h_vtx_x_L_vs_th_x_L = f4.Histo2D(ROOT::RDF::TH2DModel("h_vtx_x_L_vs_th_x_L", ";#theta_{x}^{L};x^{*,L}", 100, -600E-6, +600E-6, 100, -0.5, +0.5), "k_th_x_L", "k_vtx_x_L");
    auto h_vtx_x_R_vs_th_x_R = f4.Histo2D(ROOT::RDF::TH2DModel("h_vtx_x_R_vs_th_x_R", ";#theta_{x}^{R};x^{*,R}", 100, -600E-6, +600E-6, 100, -0.5, +0.5), "k_th_x_R", "k_vtx_x_R");
    auto h_vtx_y_L_vs_th_y_L = f4.Histo2D(ROOT::RDF::TH2DModel("h_vtx_y_L_vs_th_y_L", ";#theta_{y}^{L};y^{*,L}", 100, -600E-6, +600E-6, 100, -0.5, +0.5), "k_th_y_L", "k_vtx_y_L");
    auto h_vtx_y_R_vs_th_y_R = f4.Histo2D(ROOT::RDF::TH2DModel("h_vtx_y_R_vs_th_y_R", ";#theta_{y}^{R};y^{*,R}", 100, -600E-6, +600E-6, 100, -0.5, +0.5), "k_th_y_R", "k_vtx_y_R");

    // Line 1236:
    //           (k.vtx_x_R - k.vtx_x_L)
    //           (k.vtx_y_R - k.vtx_y_L)
    auto h_vtx_x_diffLR = f4.Histo1D(ROOT::RDF::TH1DModel("h_vtx_x_diffLR", ";x^{*,R} - x^{*,L}", 100, -0.5, +0.5), "k_vtx_x_diffLR");
    auto h_vtx_y_diffLR = f4.Histo1D(ROOT::RDF::TH1DModel("h_vtx_y_diffLR", ";y^{*,R} - y^{*,L}", 100, -0.5, +0.5), "k_vtx_y_diffLR");

    // Line 1239:
    //           (k.th_x, k.vtx_x_R - k.vtx_x_L)
    //           (k.th_y, k.vtx_y_R - k.vtx_y_L)
    auto h_vtx_x_diffLR_vs_th_x = f4.Histo1D(ROOT::RDF::TH1DModel("h_vtx_x_diffLR", ";x^{*,R} - x^{*,L}", 100, -0.5, +0.5), "k_th_x", "k_vtx_x_diffLR");
    auto h_vtx_y_diffLR_vs_th_y = f4.Histo1D(ROOT::RDF::TH1DModel("h_vtx_y_diffLR", ";y^{*,R} - y^{*,L}", 100, -0.5, +0.5), "k_th_y", "k_vtx_y_diffLR");

    // Line 1245:
    //           (k.vtx_x_R, k.vtx_x_R - k.vtx_x_L)
    //           (k.vtx_y_R, k.vtx_y_R - k.vtx_y_L)
    auto h_vtx_x_diffLR_vs_vtx_x_R = f4.Histo2D(ROOT::RDF::TH2DModel("h_vtx_x_diffLR_vs_vtx_x_R", ";x^{*,R};x^{*,R} - x^{*,L}", 100, -0.5, +0.5, 100, -0.5, +0.5), "k_vtx_x_R", "k_vtx_y_diffLR");
    auto h_vtx_y_diffLR_vs_vtx_y_R = f4.Histo2D(ROOT::RDF::TH2DModel("h_vtx_y_diffLR_vs_vtx_y_R", ";y^{*,R};y^{*,R} - y^{*,L}", 100, -0.5, +0.5, 100, -0.5, +0.5), "k_vtx_y_R", "k_vtx_y_diffLR");

    auto r7 = f4.Define("correction", "CalculateAcceptanceCorrectionsRDF(kinematics)")
                 .Define("corr",       "correction.corr")
                 .Define("div_corr",   "correction.div_corr")
                 .Define("one",        "One()");

   Binning* bis;
   for(int bi = 0; bi < binnings.size() ; bi++){
     bis = &binning_setup[bi];
     bh_t_Nev_before[bi] = r7.Histo1D(ROOT::RDF::TH1DModel("h_t_Nev_before", ";|t|;events per bin", bis->N_bins, bis->bin_edges), "k_t", "one");
     bh_t_before[bi] = r7.Histo1D(ROOT::RDF::TH1DModel("h_t_before", ";|t|", bis->N_bins, bis->bin_edges), "k_t", "one");
   }

   // Line 1412
   auto h_th_y_vs_th_x_before = r7.Histo2D(ROOT::RDF::TH2DModel("h_th_y_vs_th_x_before", ";#theta_{x};#theta_{y}", 150, -300E-6, +300E-6, 150, -150E-6, +150E-6), "k_th_x", "k_th_y", "one");

   auto f5 = r7.Filter("! correction.skip", "acceptance correction");

   for(int bi = 0; bi < binnings.size(); bi++){
     bis = &binning_setup[bi];
     bh_t_Nev_after_no_corr[bi] = f5.Histo1D(ROOT::RDF::TH1DModel("h_t_Nev_after_no_corr", ";|t|;events per bin", bis->N_bins, bis->bin_edges), "k_t", "one");
     bh_t_after_no_corr[bi] = f5.Histo1D(ROOT::RDF::TH1DModel("h_t_after_no_corr", ";|t|", bis->N_bins, bis->bin_edges), "k_t", "one");
     bh_t_after[bi] = f5.Histo1D(ROOT::RDF::TH1DModel("h_t_after", ";|t|", bis->N_bins, bis->bin_edges), "k_t", "corr");
   }

   // Line 1435
   auto h_th_y_vs_th_x_after = f5.Histo2D(ROOT::RDF::TH2DModel("h_th_y_vs_th_x_after", ";#theta_{x};#theta_{y}", 150, -300E-6, +300E-6, 150, -150E-6, +150E-6), "k_th_x", "k_th_y", "div_corr");

   // Line 1435
   auto h_th_vs_phi_after = f5.Histo2D(ROOT::RDF::TH2DModel("h_th_vs_phi_after", ";#phi;#theta", 50, -M_PI, +M_PI, 50, 150E-6, 550E-6), "k_phi", "k_th", "div_corr");

   // Line 1441
   // apply normalization
   auto bh_t_normalized_ob_1_30_02 = f5.Define("corr_norm", "corr * normalization")
                                  .Histo1D(ROOT::RDF::TH1DModel("h_t_normalized", ";|t|",128, 0., 4.), "k_t", "corr_norm");

   // Line 1445
   auto h_th_y_vs_th_x_normalized = f5.Define("div_corr_norm", "correction.div_corr * normalization")
                                 .Histo2D(ROOT::RDF::TH2DModel("h_th_y_vs_th_x_normalized", ";#theta_{x};#theta_{y}", 150, -600E-6, +600E-6, 150, -600E-6, +600E-6), "k_th_x", "k_th_y", "div_corr_norm");

   // Trigger event
   //h_y_L_1_F_vs_x_L_1_F_al_nosel.GetValue();
   auto c = f5.Count().GetValue();
   std::cout << c << std::endl;
}
