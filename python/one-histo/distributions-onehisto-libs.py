import ROOT
import sys
import os.path

RDF = ROOT.ROOT.RDataFrame

# Load C++ headers with definitions
ROOT.gInterpreter.Declare('#include "libs/distributions_lib.h"')

# Load share library with implementations
ROOT.gSystem.Load("libs/libdistributions.so")

# FIXME: Do not manage to get it from ROOT
M_PI = 3.14159265358979323846264338328      # Pi

if len(sys.argv) < 2:
    print('Usage: python distributions.py input_filename')
    sys.exit(1)  # no input file specified

if len(sys.argv) == 3:
    nthreads = int(sys.argv[2])
    if nthreads != 0:
        ROOT.ROOT.EnableImplicitMT(nthreads)

# Read input files
fname = sys.argv[1]

# Get input (temporal options)
treename          = "distilled"
selected_diagonal = "d45b_56t"
prefix            = ""
outputDir         = "."
input_file        = prefix + fname # Created with distill.py

if not os.path.isfile(input_file):
    print('File does not exists: %s' % input_file)
    sys.exit(1)

# Line 112
# init diagonal settings
ROOT.Init("45b_56t");

# Line 117
# default parameters
detailsLevel         = 0 	    # 0: no details, 1: some details, >= 2 all details
overrideCutSelection = False	# whether the default cut selection should be overriden by the command-line selection
cutSelectionString   = None
outputDir            = "."
inputDir             = "."
input_n_si           = 4.0
time_group_divisor   = 0
time_group_remainder = 0
event_group_divisor  = 0
event_group_index    = 0
evIdxStep            = 1
maxTaggedEvents      = 0	   # 0 means no maximum

# Line 131
# parse command line arguments, starting from index 2

# Line 226 - Print parameters
print("* detailsLevel = %s" % detailsLevel)
print("* outputDir = %s" % outputDir)
print("* inputDir = %s" % inputDir)
print("* input n_si = %s" % input_n_si)
print("* time_group_divisor = %s" % time_group_divisor)
print("* time_group_remainder = %s" % time_group_remainder)
print("* event_group_divisor = %s" % event_group_divisor)
print("* event_group_index = %s" % event_group_index)
print("* evIdxStep = %s" % evIdxStep)
print("* maxTaggedEvents = %s" % maxTaggedEvents)

# Line 237
# select cuts
ROOT.anal.BuildCuts()
ROOT.anal.n_si = input_n_si

# Line 241 - Let's start assuming no overrideCutSelection
# TODO if (overrideCutSelection)

# Line 272
# print info
print("\n");
print("------------------------------ environment ------------------------------\n");
ROOT.env.Print();
print("\n");
print("------------------------------- analysis --------------------------------\n");
ROOT.anal.Print();
print("\n");

# Line 281
# alignment init
for i,_ in  enumerate(ROOT.alignmentSources):
	print("\n---------- alignment source %s ----------\n" % i);
	ROOT.alignmentSources[i].Init();

print("\n\n");

# Line 289
# binnings
binnings = ROOT.vector('string')()
binnings.push_back("ub");
binnings.push_back("ob-1-10-0.2");
binnings.push_back("ob-1-30-0.2");


#########################################################
######   READ INPUT FILE, INITIALIZE RDATAFRAME    ######
#########################################################

# Line 286 - 325
#     equivalent to: get input
#                    init input data
#                    get input data
# Read all branches
rdf = RDF(treename, input_file)

# Line 327
# get time-dependent corrections
corrg_pileup = None
if ROOT.anal.use_pileup_efficiency_fits:
	path = inputDir + "/pileup_fit_combined.root"
	puF = ROOT.TFile.Open(path)
	if not os.path.exists(puF):
		print("ERROR: pile-up correction file `%s' cannot be opened.\n" % path);
	if diagonal == "d45b_56t":
		#corrg_pileup = (TGraph *) puF.Get("45b_56t/dgn");
		corrg_pileup = puF.Get("45b_56t/dgn")
	if diagonal == "d45t_56b":
		#corrg_pileup = (TGraph *) puF.Get("45b_56t/dgn");
		corrg_pileup = puF.Get("45t_56b/dgn")

# Line 358
# get th_y* dependent efficiency correction
f_3outof4_efficiency_L_F = None;
f_3outof4_efficiency_L_N = None;
f_3outof4_efficiency_R_N = None;
f_3outof4_efficiency_R_F = None;

if ROOT.anal.use_3outof4_efficiency_fits:
	path = inputDir + "/eff3outof4_details_fit_old.root"
	effFile = ROOT.TFile.Open(path)
	if (os.path.exists(effFile)):
		print("ERROR: 3-out-of-4 efficiency file `%s' cannot be opened.\n" % path);

	diagonal = selected_diagonal;
	f_3outof4_efficiency_L_F = effFile.Get(diagonal + "/L_F/fit");
	f_3outof4_efficiency_L_N = effFile.Get(diagonal + "/L_N/fit");
	f_3outof4_efficiency_R_N = effFile.Get(diagonal + "/R_N/fit");
	f_3outof4_efficiency_R_F = effFile.Get(diagonal + "/R_F/fit");

	print("\n>> using 3-out-of-4 fits: %s, %s, %s, %s\n" %
		(f_3outof4_efficiency_L_F, f_3outof4_efficiency_L_N,
         f_3outof4_efficiency_R_N, f_3outof4_efficiency_R_F))

# TODO Line 380 (not needed AFAIK)
# get unsmearing correction

# Line 394
# book metadata histograms
ROOT.timestamp_bins = ROOT.timestamp_max - ROOT.timestamp_min + 1.;

# Long TODO
# Lines 397 - 779
# THistograms, TGraphs and TProfiles declaration

# Create bh_t_* hists
# FIXME Define proper binnings

bh_t_Nev_before = dict()
bh_t_Nev_after_no_corr = dict()
bh_t_before = dict()
bh_t_after = dict()
bh_t_after_no_corr = dict()
bp_t_phi_corr = dict()
bp_t_full_corr = dict()

binning_setup = dict()
for b in binnings:
    binning = ROOT.BuildBinningRDF(ROOT.anal, b)
    binning_setup[b] = binning

#Line 780
# zero counters
n_ev_full = 0;
n_ev_cut = dict();
for ci in range(ROOT.anal.N_cuts):
	n_ev_cut[ci] = 0

th_min = 1E100;
th_y_L_min = +1E100; th_y_R_min = +1E100

N_anal=0; N_anal_zeroBias=0;
N_zeroBias_el=0; N_zeroBias_el_RP_trig=0;
N_4outof4=0; N_el=0;
N_el_T2trig=0; N_4outof4_T2trig=0;
N_el_raw=0;

# TODO
# Line 795
# map<unsigned int, pair<unsigned int, unsigned int> > runTimestampBoundaries;


#########################################################
###### FILTER, BUILD HISTOGRAMS - START EVENT LOOP ######
#########################################################


f1 = rdf.Filter("! SkipTime( timestamp )", 'check time - selected')

# Line 822
if time_group_divisor != 0:
    f1 = f1.Filter("! SkipTimeInterval( timestamp, %s, %s )".format(time_group_divisor, time_group_remainder),
                    'time interval')

# Diagonal cut (L831)
f2 = f1.Filter("v_L_2_F && v_L_2_N && v_R_2_F && v_R_2_N", 'allDiagonalRPs')

# Not cut for this filter in original code
f_zerobias = f2.Filter("! ((trigger_bits & 512) != 0)", 'zero_bias_event')

xs = ["x_L_1_F", "x_L_2_N", "x_L_2_F", "x_R_1_F", "x_R_2_N", "x_R_2_F"]
ys = ["y_L_1_F", "y_L_2_N", "y_L_2_F", "y_R_1_F", "y_R_2_N", "y_R_2_F"]

# Apply fine alignment (L 852)
r2 = f2.Define("h_al", "ApplyFineAlignment( timestamp ,{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {})".format(*(xs+ys) )) \
       .Define("h_al_x_L_1_F", "h_al.L_1_F.x") \
       .Define("h_al_x_L_2_N", "h_al.L_2_N.x") \
       .Define("h_al_x_L_2_F", "h_al.L_2_F.x") \
       .Define("h_al_y_L_1_F", "h_al.L_1_F.y") \
       .Define("h_al_y_L_2_N", "h_al.L_2_N.y") \
       .Define("h_al_y_L_2_F", "h_al.L_2_F.y") \
       .Define("h_al_x_R_1_F", "h_al.R_1_F.x") \
       .Define("h_al_x_R_2_N", "h_al.R_2_N.x") \
       .Define("h_al_x_R_2_F", "h_al.R_2_F.x") \
       .Define("h_al_y_R_1_F", "h_al.R_1_F.y") \
       .Define("h_al_y_R_2_N", "h_al.R_2_N.y") \
       .Define("h_al_y_R_2_F", "h_al.R_2_F.y") \

# fill pre-selection histograms (Line 860 - 866)
#
#al_nosel_models = map(ROOT.ROOT.RDF.TH2DModel, [
#    ("h_y_L_1_F_vs_x_L_1_F_al_nosel", ";x^{L,1,F};y^{L,1,F}", 150, -15., 15., 300, -30., +30.),
#    ("h_y_L_2_N_vs_x_L_2_N_al_nosel", ";x^{L,2,N};y^{L,2,N}", 150, -15., 15., 300, -30., +30.),
#    ("h_y_L_2_F_vs_x_L_2_F_al_nosel", ";x^{L,2,F};y^{L,2,F}", 150, -15., 15., 300, -30., +30.),
#    ("h_y_R_1_F_vs_x_R_1_F_al_nosel", ";x^{R,1,F};y^{R,1,F}", 150, -15., 15., 300, -30., +30.),
#    ("h_y_R_2_N_vs_x_R_2_N_al_nosel", ";x^{R,2,N};y^{R,2,N}", 150, -15., 15., 300, -30., +30.),
#    ("h_y_R_2_F_vs_x_R_2_F_al_nosel", ";x^{R,2,F};y^{R,2,F}", 150, -15., 15., 300, -30., +30.)
#])

al_nosel_models0 = ROOT.ROOT.RDF.TH2DModel(("h_y_L_1_F_vs_x_L_1_F_al_nosel", ";x^{L,1,F};y^{L,1,F}", 150, -15., 15., 300, -30., +30.)) 
h_y_L_1_F_vs_x_L_1_F_al_nosel = r2.Histo2D(al_nosel_models0, "h_al_x_L_1_F", "h_al_y_L_1_F")
#h_y_L_2_N_vs_x_L_2_N_al_nosel = r2.Histo2D(al_nosel_models[1], "h_al_x_L_2_N", "h_al_y_L_2_N")
#h_y_L_2_F_vs_x_L_2_F_al_nosel = r2.Histo2D(al_nosel_models[2], "h_al_x_L_2_F", "h_al_y_L_2_F")
#h_y_R_1_F_vs_x_R_1_F_al_nosel = r2.Histo2D(al_nosel_models[3], "h_al_x_R_1_F", "h_al_y_R_1_F")
#h_y_R_2_N_vs_x_R_2_N_al_nosel = r2.Histo2D(al_nosel_models[4], "h_al_x_R_2_N", "h_al_y_R_2_N")
#h_y_R_2_F_vs_x_R_2_F_al_nosel = r2.Histo2D(al_nosel_models[5], "h_al_x_R_2_F", "h_al_y_R_2_F")

# run reconstruction (Line 876)
### kinematics struct
r3 = r2.Define("kinematics", 'DoReconstruction( h_al )')

r4 = r3.Define("k_th_x_R",        "kinematics.th_x_R") \
       .Define("k_th_y_R",        "kinematics.th_y_R") \
       .Define("k_th_x_L",        "kinematics.th_x_L") \
       .Define("k_th_y_L",        "kinematics.th_y_L") \
       .Define("k_th_x",          "kinematics.th_x") \
       .Define("k_th_y",          "kinematics.th_y") \
       .Define("minus_k_th_y",    "- kinematics.th_y") \
       .Define("k_vtx_x",         "kinematics.vtx_x") \
       .Define("k_vtx_x_L",       "kinematics.vtx_x_L") \
       .Define("k_vtx_x_R",       "kinematics.vtx_x_R") \
       .Define("k_vtx_y",         "kinematics.vtx_y")   \
       .Define("k_vtx_y_L",       "kinematics.vtx_y_L") \
       .Define("k_vtx_y_R",       "kinematics.vtx_y_R") \
       .Define("k_th_y_L_F",      "kinematics.th_y_L_F") \
       .Define("k_th_y_L_N",      "kinematics.th_y_L_N") \
       .Define("k_th_y_R_F",      "kinematics.th_y_R_F") \
       .Define("k_th_y_R_N",      "kinematics.th_y_R_N") \
       .Define("k_th_x_diffLR",   "kinematics.th_x_R - kinematics.th_x_L") \
       .Define("k_th_y_diffLR",   "kinematics.th_y_R - kinematics.th_y_L") \
       .Define("k_th_x_diffLF",   "kinematics.th_x_L - kinematics.th_x") \
       .Define("k_th_x_diffRF",   "kinematics.th_x_R - kinematics.th_x") \
       .Define("k_th_y_L_diffNF", "kinematics.th_y_L_F - kinematics.th_y_L_N") \
       .Define("k_th_y_R_diffNF", "kinematics.th_y_R_F - kinematics.th_y_R_N") \
       .Define("k_vtx_x_diffLR",  "kinematics.vtx_x_R - kinematics.vtx_x_L") \
       .Define("k_vtx_y_diffLR",  "kinematics.vtx_y_R - kinematics.vtx_y_L") \
       .Define("k_t",             "kinematics.t")                            \
       .Define("k_th",            "kinematics.th")                           \
       .Define("k_phi",           "kinematics.phi")

# cut evaluation
r5 = r4.Define("cutdata", "EvaluateCutsRDF( h_al, kinematics )")

# Elastic cut
f4 = r5.Filter("cutdata.select", "elastic cut")

# Define normalization and norm_corr colums
r6 = r5.Define("norm_corr",     "getNorm_corr( timestamp )" ) \
       .Define("normalization", "getNormalization( norm_corr )")

# Line 1401
# calculate acceptance divergence correction
r7 = f4.Define("correction", "CalculateAcceptanceCorrectionsRDF( kinematics )") \
       .Define("corr",       "correction.corr") \
       .Define("div_corr",   "correction.div_corr") \
       .Define("one",        "One()")

# Line 1414
# Filter skip
f5 = r7.Filter("! correction.skip", "acceptance correction")

###############################
###### END OF EVENT LOOP ######
###############################

# Trigger event
h_y_L_1_F_vs_x_L_1_F_al_nosel.GetValue();
