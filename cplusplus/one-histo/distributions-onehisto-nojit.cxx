#include <ROOT/RDataFrame.hxx>

#include <TInterpreter.h>

#include "../common_definitions.h"
#include "../parameters_global.h"
#include "../common_algorithms.h"
#include "../parameters.h"
#include "../common.h"

#include <iostream>
#include <cstdlib>

#define GETMEMBER( type, member ) [](type &st){ return st.member; }
#define GETEXPR( type, m1, symbol, m2 ) [](type &st){ return st.m1 symbol st.m2 ;}
#define GETINVERSE( type, member) [](type &st){ return - st.member; }

using RDF = ROOT::RDataFrame;
using namespace std;

// Get input (temporal options)
auto treename          = "distilled";
auto selected_diagonal = "d45b_56t";
auto prefix            = "";
auto outputDir         = ".";

int main(int argc, char **argv)
{
  // Enable implicit parallelism
  if(argc > 2 && atoi(argv[2]) != 0)
    ROOT::EnableImplicitMT(atoi(argv[2]));

  gInterpreter->Declare(R"cpp(
       #include "../common_definitions.h"
       #include "../parameters_global.h"
       #include "../common_algorithms.h"
       #include "../parameters.h"
       #include "../common.h"
  )cpp");
  // Read input files
  auto fname = argv[1];
  auto input_file  = fname ; // Created with distill.py

  Init("45b_56t");

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
  anal.BuildCuts();
  anal.n_si = input_n_si;

  // print info
  printf("\n");
  printf("------------------------------ environment ------------------------------\n");
  env.Print();
  printf("\n");
  printf("------------------------------- analysis --------------------------------\n");
  anal.Print();
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

  for (unsigned int bi = 0; bi < binnings.size(); ++bi)
  {
    Binning b;
    BuildBinningRDF(anal, binnings[bi], b);
    binning_setup[bi] = b;
  }
  // zero counters
     unsigned long n_ev_full = 0;
     map<unsigned int, unsigned long> n_ev_cut;
     for (unsigned int ci = 1; ci <= anal.N_cuts; ++ci)
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


  auto SkipTimeBis = [](unsigned int &t){
    return ! SkipTime(t);
  };

  auto f1 = rdf.Filter( SkipTimeBis, {"timestamp"}, "check time - selected");

  // Diagonal cut (L831)
  auto allDiagonalRPs = [](unsigned int &v_L_2_F,unsigned int &v_L_2_N,unsigned int &v_R_2_F,unsigned int &v_R_2_N){
    return v_L_2_F && v_L_2_N && v_R_2_F && v_R_2_N;
  };

  auto f2 = f1.Filter(allDiagonalRPs, {"v_L_2_F", "v_L_2_N", "v_R_2_F", "v_R_2_N"}, "allDiagonalRPs");

  // Not cut for this filter in original code
  auto isZeroBiasEvent = [](unsigned int &bits){
      return ! ((bits & 512) != 0);
  };

  auto f_zerobias = f2.Filter(isZeroBiasEvent, {"trigger_bits"}, "zero_bias_event");

  auto r2 = f2.Define("h_al", ApplyFineAlignment,  {"timestamp" , "x_L_1_F", "x_L_2_N", "x_L_2_F", "x_R_1_F", "x_R_2_N", "x_R_2_F", "y_L_1_F", "y_L_2_N", "y_L_2_F", "y_R_1_F", "y_R_2_N", "y_R_2_F"})
              .Define("h_al_x_L_1_F", GETMEMBER(HitData, L_1_F.x), {"h_al"})
              .Define("h_al_x_L_2_N", GETMEMBER(HitData, L_2_N.x), {"h_al"})
              .Define("h_al_x_L_2_F", GETMEMBER(HitData, L_2_F.x), {"h_al"})
              .Define("h_al_y_L_1_F", GETMEMBER(HitData, L_1_F.y), {"h_al"})
              .Define("h_al_y_L_2_N", GETMEMBER(HitData, L_2_N.y), {"h_al"})
              .Define("h_al_y_L_2_F", GETMEMBER(HitData, L_2_F.y), {"h_al"})
              .Define("h_al_x_R_1_F", GETMEMBER(HitData, R_1_F.x), {"h_al"})
              .Define("h_al_x_R_2_N", GETMEMBER(HitData, R_2_N.x), {"h_al"})
              .Define("h_al_x_R_2_F", GETMEMBER(HitData, R_2_F.x), {"h_al"})
              .Define("h_al_y_R_1_F", GETMEMBER(HitData, R_1_F.y), {"h_al"})
              .Define("h_al_y_R_2_N", GETMEMBER(HitData, R_2_N.y), {"h_al"})
              .Define("h_al_y_R_2_F", GETMEMBER(HitData, R_2_F.y), {"h_al"});

//   // fill pre-selection histograms (Line 860 - 866)
   // al_nosel_models = {
   //    ("h_y_L_1_F_vs_x_L_1_F_al_nosel", ";x^{L,1,F};y^{L,1,F}", 150, -15., 15., 300, -30., +30.),
   //    ("h_y_L_2_N_vs_x_L_2_N_al_nosel", ";x^{L,2,N};y^{L,2,N}", 150, -15., 15., 300, -30., +30.),
   //    ("h_y_L_2_F_vs_x_L_2_F_al_nosel", ";x^{L,2,F};y^{L,2,F}", 150, -15., 15., 300, -30., +30.),
   //    ("h_y_R_1_F_vs_x_R_1_F_al_nosel", ";x^{R,1,F};y^{R,1,F}", 150, -15., 15., 300, -30., +30.),
   //    ("h_y_R_2_N_vs_x_R_2_N_al_nosel", ";x^{R,2,N};y^{R,2,N}", 150, -15., 15., 300, -30., +30.),
   //    ("h_y_R_2_F_vs_x_R_2_F_al_nosel", ";x^{R,2,F};y^{R,2,F}", 150, -15., 15., 300, -30., +30.)
   //};

   auto al_nosel_models0 = ROOT::RDF::TH2DModel("h_y_L_1_F_vs_x_L_1_F_al_nosel", ";x^{L,1,F};y^{L,1,F}", 150, -15., 15., 300, -30., +30.);
   auto h_y_L_1_F_vs_x_L_1_F_al_nosel = r2.Histo2D(al_nosel_models0, "h_al_x_L_1_F", "h_al_y_L_1_F");
   //auto h_y_L_2_N_vs_x_L_2_N_al_nosel = r2.Histo2D(al_nosel_models[1], "h_al_x_L_2_N", "h_al_y_L_2_N");
   //auto h_y_L_2_F_vs_x_L_2_F_al_nosel = r2.Histo2D(al_nosel_models[2], "h_al_x_L_2_F", "h_al_y_L_2_F");
   //auto h_y_R_1_F_vs_x_R_1_F_al_nosel = r2.Histo2D(al_nosel_models[3], "h_al_x_R_1_F", "h_al_y_R_1_F");
   //auto h_y_R_2_N_vs_x_R_2_N_al_nosel = r2.Histo2D(al_nosel_models[4], "h_al_x_R_2_N", "h_al_y_R_2_N");
   //auto h_y_R_2_F_vs_x_R_2_F_al_nosel = r2.Histo2D(al_nosel_models[5], "h_al_x_R_2_F", "h_al_y_R_2_F");

   auto r3 = r2.Define("kinematics", DoReconstruction , { "h_al" } );

   auto r4 = r3.Define("k_th_x_R",         GETMEMBER( Kinematics, th_x_R ), {"kinematics"})
                .Define("k_th_y_R",        GETMEMBER( Kinematics, th_y_R), {"kinematics"})
                .Define("k_th_x_L",        GETMEMBER( Kinematics, th_x_L), {"kinematics"})
                .Define("k_th_y_L",        GETMEMBER( Kinematics, th_y_L), {"kinematics"})
                .Define("k_th_x",          GETMEMBER( Kinematics, th_x), {"kinematics"})
                .Define("k_th_y",          GETMEMBER( Kinematics, th_y), {"kinematics"})
                .Define("minus_k_th_y",    GETINVERSE(Kinematics, th_y), {"kinematics"})
                .Define("k_vtx_x",         GETMEMBER( Kinematics, vtx_x), {"kinematics"})
                .Define("k_vtx_x_L",       GETMEMBER( Kinematics, vtx_x_L), {"kinematics"})
                .Define("k_vtx_x_R",       GETMEMBER( Kinematics, vtx_x_R), {"kinematics"})
                .Define("k_vtx_y",         GETMEMBER( Kinematics, vtx_y), {"kinematics"})
                .Define("k_vtx_y_L",       GETMEMBER( Kinematics, vtx_y_L), {"kinematics"})
                .Define("k_vtx_y_R",       GETMEMBER( Kinematics, vtx_y_R), {"kinematics"})
                .Define("k_th_y_L_F",      GETMEMBER( Kinematics, th_y_L_F), {"kinematics"})
                .Define("k_th_y_L_N",      GETMEMBER( Kinematics, th_y_L_N), {"kinematics"})
                .Define("k_th_y_R_F",      GETMEMBER( Kinematics, th_y_R_F), {"kinematics"})
                .Define("k_th_y_R_N",      GETMEMBER( Kinematics, th_y_R_N), {"kinematics"})
                .Define("k_th_x_diffLR",   GETEXPR(   Kinematics, th_x_R, - , th_x_L), {"kinematics"})
                .Define("k_th_y_diffLR",   GETEXPR(   Kinematics, th_y_R, - , th_y_L), {"kinematics"})
                .Define("k_th_x_diffLF",   GETEXPR(   Kinematics, th_x_L, - , th_x), {"kinematics"})
                .Define("k_th_x_diffRF",   GETEXPR(   Kinematics, th_x_R, - , th_x), {"kinematics"})
                .Define("k_th_y_L_diffNF", GETEXPR(   Kinematics, th_y_L_F, - , th_y_L_N), {"kinematics"})
                .Define("k_th_y_R_diffNF", GETEXPR(   Kinematics, th_y_R_F, - , th_y_R_N), {"kinematics"})
                .Define("k_vtx_x_diffLR",  GETEXPR(   Kinematics, vtx_x_R, - , vtx_x_L), {"kinematics"})
                .Define("k_vtx_y_diffLR",  GETEXPR(   Kinematics, vtx_y_R, - , vtx_y_L), {"kinematics"})
                .Define("k_t",             GETMEMBER( Kinematics, t), {"kinematics"})
                .Define("k_th",            GETMEMBER( Kinematics, th), {"kinematics"})
                .Define("k_phi",           GETMEMBER( Kinematics, phi), {"kinematics"});

    auto r5 = r4.Define("cutdata", EvaluateCutsRDF, { "h_al", "kinematics"});

    // Elastic cut
    auto f4 = r5.Filter( [](CutData &cutdata){return cutdata.select;}, {"cutdata"}, "elastic cut");

    // Define normalization and norm_corr colums
    auto r6 = f4.Define("norm_corr",     getNorm_corr, {"timestamp"} )
                .Define("normalization", getNormalization, {"norm_corr"} );

    auto r7 = r6.Define("correction", CalculateAcceptanceCorrectionsRDF, {"kinematics"})
                 .Define("corr",      GETMEMBER(Correction, corr),  {"correction"})
                 .Define("div_corr",  GETMEMBER(Correction, div_corr), {"correction"})
                 .Define("one",       [](){return 1;});

    auto skipCorrection = [](Correction &correction){
        return ! correction.skip;
    };

    auto f5 = r7.Filter( skipCorrection, {"correction"}, "acceptance correction");

   // Trigger event
   h_y_L_1_F_vs_x_L_1_F_al_nosel.GetValue();
}
