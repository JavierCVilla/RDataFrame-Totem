import ROOT
import sys
import os.path

RDF = ROOT.ROOT.RDataFrame

# Load C++ headers
ROOT.gInterpreter.Declare('#include "common_definitions.h"')
ROOT.gInterpreter.Declare('#include "parameters_global.h"')
ROOT.gInterpreter.Declare('#include "common_algorithms.h"')
ROOT.gInterpreter.Declare('#include "parameters.h"')
ROOT.gInterpreter.Declare('#include "common.h"')

# FIXME: Do not manage to get it from ROOT
M_PI = 3.14159265358979323846264338328      # Pi

if len(sys.argv) < 2:
    print('Usage: python distributions.py input_filename')
    sys.exit(1)  # no input file specified

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

for bi in binnings:
    pass

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


# TODO
# Line 802

# Initial cuts

# remove troublesome runs
# unsigned int run = ev.run_num / 100000;
# unsigned int file = ev.run_num % 100000;

# NOT IMPLEMENTED: return always false
# if (SkipRun(run, file, true))
# 	continue;

# USED ANYWHERE?
# update timestamp run boundaries
# auto rtbit = runTimestampBoundaries.find(run);
# if (rtbit == runTimestampBoundaries.end())
# {
# 	runTimestampBoundaries.insert({run, {ev.timestamp, ev.timestamp}});
# } else {
# 	rtbit->second.first = min(rtbit->second.first, ev.timestamp);
# 	rtbit->second.second = max(rtbit->second.second, ev.timestamp);
# }

# Line 818
# check time - selected?
# if (anal.SkipTime(ev.timestamp))
# 	continue;
skipTime_code = """
bool SkipTime( unsigned long long &timestamp){
    extern Analysis anal ;
    return anal.SkipTime(timestamp);
}
"""
ROOT.gInterpreter.Declare(skipTime_code)

f1 = rdf.Filter("! SkipTime( timestamp )", 'check time - selected')

skipTime_interval_code = """
bool SkipTimeInterval( unsigned long long &timestamp, int &tgd, int &tgr ){
    double time_group_interval = 1.;	// s
    int time_group = int(timestamp / time_group_interval);
    return  ( (time_group % tgd) != tgr);
}
"""
ROOT.gInterpreter.Declare(skipTime_interval_code)

# Line 822
if time_group_divisor != 0:
    f1 = f1.Filter("! SkipTimeInterval( timestamp, %s, %s )".format(time_group_divisor, time_group_remainder),
                    'time interval')

# Diagonal cut (L831)
f2 = f1.Filter("v_L_2_F && v_L_2_N && v_R_2_F && v_R_2_N", 'allDiagonalRPs')

# FIXME: N_4outof4_T2trig is not used to produce any plot
# Line 835
# if ((ev.trigger_bits & 64) != 0)
# 	N_4outof4_T2trig++;

# Line 838
h_timestamp_dgn = f2.Histo1D("timestamp")

# Line 841
# NOT IMPLEMENTED
# select the elastic-trigger bunch(es) only
# if (SkipBunch(run, ev.bunch_num))
#
# SkipBunch will return False for whatever DS in use,
# since all of them define keepAllBunches = True within
# parameters.h

# zero bias event?
# Create named filter with number of zero bias events
f3 = f2.Filter("! ((trigger_bits & 512) != 0)", 'zero_bias_event')

xs = ["x_L_1_F", "x_L_2_N", "x_L_2_F", "x_R_1_F", "x_R_2_N", "x_R_2_F"]
ys = ["y_L_1_F", "y_L_2_N", "y_L_2_F", "y_R_1_F", "y_R_2_N", "y_R_2_F"]

# Apply fine alignment (L 852)
r2 = f3.Define("h_al", "ApplyFineAlignment( timestamp ,{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {})".format(*(xs+ys) )) \
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

al_nosel_models = [
    ("h_y_L_1_F_vs_x_L_1_F_al_nosel", ";x^{L,1,F};y^{L,1,F}", 150, -15., 15., 300, -30., +30.),
    ("h_y_L_2_N_vs_x_L_2_N_al_nosel", ";x^{L,2,N};y^{L,2,N}", 150, -15., 15., 300, -30., +30.),
    ("h_y_L_2_F_vs_x_L_2_F_al_nosel", ";x^{L,2,F};y^{L,2,F}", 150, -15., 15., 300, -30., +30.),
    ("h_y_R_1_F_vs_x_R_1_F_al_nosel", ";x^{R,1,F};y^{R,1,F}", 150, -15., 15., 300, -30., +30.),
    ("h_y_R_2_N_vs_x_R_2_N_al_nosel", ";x^{R,2,N};y^{R,2,N}", 150, -15., 15., 300, -30., +30.),
    ("h_y_R_2_F_vs_x_R_2_F_al_nosel", ";x^{R,2,F};y^{R,2,F}", 150, -15., 15., 300, -30., +30.)
]

h_y_L_1_F_vs_x_L_1_F_al_nosel = r2.Histo2D(al_nosel_models[0], "h_al_x_L_1_F", "h_al_y_L_1_F")
h_y_L_2_N_vs_x_L_2_N_al_nosel = r2.Histo2D(al_nosel_models[1], "h_al_x_L_2_N", "h_al_y_L_2_N")
h_y_L_2_F_vs_x_L_2_F_al_nosel = r2.Histo2D(al_nosel_models[2], "h_al_x_L_2_F", "h_al_y_L_2_F")
h_y_R_1_F_vs_x_R_1_F_al_nosel = r2.Histo2D(al_nosel_models[3], "h_al_x_R_1_F", "h_al_y_R_1_F")
h_y_R_2_N_vs_x_R_2_N_al_nosel = r2.Histo2D(al_nosel_models[4], "h_al_x_R_2_N", "h_al_y_R_2_N")
h_y_R_2_F_vs_x_R_2_F_al_nosel = r2.Histo2D(al_nosel_models[5], "h_al_x_R_2_F", "h_al_y_R_2_F")

# TODO (not needed so far AFAIK)
#if (detailsLevel >= 2)
#{
#    h_timestamp_B0->Fill(ev.timestamp);
#    g_run_vs_timestamp->SetPoint(g_run_vs_timestamp->GetN(), ev.timestamp, ev.run_num);
#    //g_ev_num_vs_timestamp->SetPoint(g_ev_num_vs_timestamp->GetN(), ev.timestamp, ev.event_num);
#    //g_tr_num_vs_timestamp->SetPoint(g_tr_num_vs_timestamp->GetN(), ev.timestamp, ev.trigger_num);
#}

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
       .Define("k_t",             "kinematics.t")

# Line 906
# cut evaluation
r5 = r4.Define("cutdata", "EvaluateCutsRDF( h_al, kinematics )")

# Line 919
# fill background distributions
# Not IMPLEMENTED, not used ANYWHERE

## TODO fill no-cut histograms (L957)
# for (unsigned int ci = 1; ci <= anal.N_cuts; ++ci)
# {
# 	//h2_cq_full[ci]->Fill(ccb[ci]*cqa[ci] - cca[ci]*cqb[ci], cca[ci]*cqa[ci] + ccb[ci]*cqb[ci] + ccc[ci]);
# 	h2_cq_full[ci]->Fill(cd.cqa[ci], cd.cqb[ci]);
# }

# Line 979
# Elastic cut
f4 = r5.Filter("cutdata.select", "elastic cut")

# TODO (L985)
# if (maxTaggedEvents > 0 && N_el > maxTaggedEvents)
# 		break;

# TODO (L988)
# g_selected_bunch_num_vs_timestamp->SetPoint(g_selected_bunch_num_vs_timestamp->GetN(), ev.timestamp, ev.bunch_num);

# TODO Line 1010
# if event_group_divisor > 0:
# 	int event_group = (N_el_raw-1) / event_group_divisor;
# 	if (event_group < event_group_index)
# 		continue;
# 	if (event_group > event_group_index)
# 		break;
# }

## TODO : Not used AFAIK
## use_3outof4_efficiency_fits

# Line 1039
getNorm_corr_code = """
double getNorm_corr( unsigned long long &timestamp ){
    extern Analysis anal;

    // determine normalization factors (luminosity + corrections)
	double inefficiency_3outof4 = anal.inefficiency_3outof4;
    double inefficiency_shower_near = anal.inefficiency_shower_near;
    double inefficiency_pile_up = anal.inefficiency_pile_up;
    double inefficiency_trigger = anal.inefficiency_trigger;

    double norm_corr =
		1./(1. - (inefficiency_3outof4 + inefficiency_shower_near))
		* 1./(1. - inefficiency_pile_up)
		* 1./(1. - inefficiency_trigger);

    return norm_corr;
}
"""
ROOT.gInterpreter.Declare(getNorm_corr_code)

# Line 1048
getNormalization_code = """
double getNormalization( double &norm_corr ){
    extern Analysis anal;

    double normalization = anal.bckg_corr * norm_corr / anal.L_int;

    return normalization;
}
"""
ROOT.gInterpreter.Declare(getNormalization_code)

# Define normalization and norm_corr colums
r6 = r5.Define("norm_corr",     "getNorm_corr( timestamp )" ) \
       .Define("normalization", "getNormalization( norm_corr )")

# Line 1044 TODO (TProiles)
# p_norm_corr->Fill(ev.timestamp, norm_corr);
# p_3outof4_corr->Fill(k.th_y, inefficiency_3outof4);

# TODO Skipped
# Line 1050 - 1087
# // data for alignment
# // (SHOULD use hit positions WITHOUT alignment corrections, i.e. ev.h)
# signed int period = int((ev.timestamp - anal.alignment_t0) / anal.alignment_ts);


# Line 1089
# Fill raw histograms
h_timestamp_sel = f4.Histo1D("timestamp");


# Line 1110
# fill histograms
noal_sel_models = [
    ("h_y_L_1_F_vs_x_L_1_F_noal_sel", ";x^{L,1,F};y^{L,1,F}", 100, -3., +3., 300, -30., +30.),
    ("h_y_L_2_N_vs_x_L_2_N_noal_sel", ";x^{L,2,N};y^{L,2,N}", 100, -3., +3., 300, -30., +30.),
    ("h_y_L_2_F_vs_x_L_2_F_noal_sel", ";x^{L,2,F};y^{L,2,F}", 100, -3., +3., 300, -30., +30.),
    ("h_y_R_1_F_vs_x_R_1_F_noal_sel", ";x^{R,1,F};y^{R,1,F}", 100, -3., +3., 300, -30., +30.),
    ("h_y_R_2_N_vs_x_R_2_N_noal_sel", ";x^{R,2,N};y^{R,2,N}", 100, -3., +3., 300, -30., +30.),
    ("h_y_R_2_F_vs_x_R_2_F_noal_sel", ";x^{R,2,F};y^{R,2,F}", 100, -3., +3., 300, -30., +30.)
]

h_y_L_1_F_vs_x_L_1_F_noal_sel = f4.Histo2D(noal_sel_models[0], "x_L_1_F", "y_L_1_F")
h_y_L_2_N_vs_x_L_2_N_noal_sel = f4.Histo2D(noal_sel_models[1], "x_L_2_N", "y_L_2_N")
h_y_L_2_F_vs_x_L_2_F_noal_sel = f4.Histo2D(noal_sel_models[2], "x_L_2_F", "y_L_2_F")
h_y_R_1_F_vs_x_R_1_F_noal_sel = f4.Histo2D(noal_sel_models[3], "x_R_1_F", "y_R_1_F")
h_y_R_2_N_vs_x_R_2_N_noal_sel = f4.Histo2D(noal_sel_models[4], "x_R_2_N", "y_R_2_N")
h_y_R_2_F_vs_x_R_2_F_noal_sel = f4.Histo2D(noal_sel_models[5], "x_R_2_F", "y_R_2_F")


# Line 1117
al_sel_models = [
    ("h_y_L_1_F_vs_x_L_1_F_al_sel", ";x^{L,1,F};y^{L,1,F}", 100, -3., +3., 300, -30., +30.),
    ("h_y_L_2_N_vs_x_L_2_N_al_sel", ";x^{L,2,N};y^{L,2,N}", 100, -3., +3., 300, -30., +30.),
    ("h_y_L_2_F_vs_x_L_2_F_al_sel", ";x^{L,2,F};y^{L,2,F}", 100, -3., +3., 300, -30., +30.),
    ("h_y_R_1_F_vs_x_R_1_F_al_sel", ";x^{R,1,F};y^{R,1,F}", 100, -3., +3., 300, -30., +30.),
    ("h_y_R_2_N_vs_x_R_2_N_al_sel", ";x^{R,2,N};y^{R,2,N}", 100, -3., +3., 300, -30., +30.),
    ("h_y_R_2_F_vs_x_R_2_F_al_sel", ";x^{R,2,F};y^{R,2,F}", 100, -3., +3., 300, -30., +30.)
]

h_y_L_1_F_vs_x_L_1_F_al_sel = f4.Histo2D(al_sel_models[0], "x_L_1_F", "y_L_1_F")
h_y_L_2_N_vs_x_L_2_N_al_sel = f4.Histo2D(al_sel_models[1], "x_L_2_N", "y_L_2_N")
h_y_L_2_F_vs_x_L_2_F_al_sel = f4.Histo2D(al_sel_models[2], "x_L_2_F", "y_L_2_F")
h_y_R_1_F_vs_x_R_1_F_al_sel = f4.Histo2D(al_sel_models[3], "x_R_1_F", "y_R_1_F")
h_y_R_2_N_vs_x_R_2_N_al_sel = f4.Histo2D(al_sel_models[4], "x_R_2_N", "y_R_2_N")
h_y_R_2_F_vs_x_R_2_F_al_sel = f4.Histo2D(al_sel_models[5], "x_R_2_F", "y_R_2_F")

# Line 1157 (k.th_x_R - k.th_x_L)
#           (k.th_y_R - k.th_y_L)
th_x_diffLR = f4.Histo1D("k_th_x_diffLR")
th_y_diffLR = f4.Histo1D("k_th_y_diffLR")

# Line 1160 (k.th_x_L - k.th_x)
#           (k.th_x_R - k.th_x)
th_x_diffLF = f4.Histo1D("k_th_x_diffLF")
th_x_diffRF = f4.Histo1D("k_th_x_diffRF")

# Line 1163 (k.th_x, k.th_x_R - k.th_x_L)
#           (k.th_y, k.th_y_R - k.th_y_L)
#           (k.vtx_x, k.th_x_R - k.th_x_L)
h_th_x_diffLR_vs_th_x  = f4.Histo1D("k_th_x", "k_th_x_diffLR")
h_th_y_diffLR_vs_th_y  = f4.Histo1D("k_th_y", "k_th_y_diffLR")
h_th_x_diffLR_vs_vtx_x = f4.Histo1D("k_vtx_x", "k_th_x_diffLR")

# TODO:
# Line 1168 TProfile
#p_th_x_diffLR_vs_th_x = f4.Histo1D("k_th_x", "k_th_x_diffLR")
#p_th_y_diffLR_vs_th_y = f4.Histo1D("k_th_y", "k_th_y_diffLR")
#p_th_y_L_diffNF_vs_th_y_L = f4.Histo1D("k_th_y_L", "k_th_y_L_diffNF")
#p_th_y_R_diffNF_vs_th_y_R = f4.Histo1D("k_th_y_R", "k_th_y_R_diffNF")

# TODO:
# Line 1172 TProfile
#p_th_x_diffLR_vs_vtx_x = f4.Histo1D("k_vtx_x", "k_th_x_diffRF")

# TODO:
# if (safe)
# {
# 	th_y_diffLR_safe->Fill(k.th_y_R - k.th_y_L);
# 	th_x_diffLR_safe->Fill(k.th_x_R - k.th_x_L);
# }

# TODO:
# Line 1184 TProfile
#p_th_x_vs_th_y = f4.Histo1D("k_th_y", "k_th_x")
#p_th_x_L_vs_th_y_L = f4.Histo1D("k_th_y_L", "k_th_x_L")
#p_th_x_R_vs_th_y_R = f4.Histo1D("k_th_y_R", "k_th_x_R")

# Line 1188 (k.th_x_L, k.th_y_L)
#           (k.th_x_R, k.th_y_R)
#           (k.th_x, k.th_y)
h_th_y_L_vs_th_x_L = f4.Histo1D("k_th_x_L", "k_th_y_L")
h_th_y_R_vs_th_x_R = f4.Histo1D("k_th_x_R", "k_th_y_R")
h_th_y_vs_th_x     = f4.Histo1D("k_th_x", "k_th_y")

# TODO: TGraph
# Line 1192
# if (detailsLevel >= 1)
# {
# 	g_th_y_L_vs_th_x_L->SetPoint(g_th_y_L_vs_th_x_L->GetN(), k.th_x_L, k.th_y_L);
# 	g_th_y_R_vs_th_x_R->SetPoint(g_th_y_R_vs_th_x_R->GetN(), k.th_x_R, k.th_y_R);
# 	g_th_y_vs_th_x->SetPoint(g_th_y_vs_th_x->GetN(), k.th_x, k.th_y);
# }

# Line 1199 (k.th_y_R, k.th_y_L)
h_th_y_L_vs_th_y_R = f4.Histo1D("k_th_y_R", "k_th_y_L")

# TODO TGraph
# if (detailsLevel > 2)
# 	g_th_y_L_vs_th_y_R->SetPoint(g_th_y_L_vs_th_y_R->GetN(), k.th_y_R, k.th_y_L);

# Line 1203: (k.th_x)
#            (k.th_y)
h_th_x     = f4.Histo1D("k_th_x")
h_th_y     = f4.Histo1D("k_th_y")

# Line 1205: (-k.th_y)
h_th_y_flipped = f4.Histo1D("minus_k_th_y")

# Line 1207: (k.th_x_L)
#            (k.th_x_R)
h_th_x_L   = f4.Histo1D("k_th_x_L")
h_th_x_R   = f4.Histo1D("k_th_x_R")

# Line 1210: (k.th_y_L)
#            (k.th_y_R)
h_th_y_L   = f4.Histo1D("k_th_y_L")
h_th_y_R   = f4.Histo1D("k_th_y_R")

# Line 1213: (k.th_y_L_F)
#            (k.th_y_L_N)
#            (k.th_y_R_N)
#            (k.th_y_R_F)
h_th_y_L_F = f4.Histo1D("k_th_y_L_F")
h_th_y_L_N = f4.Histo1D("k_th_y_L_N")
h_th_y_R_N = f4.Histo1D("k_th_y_R_N")
h_th_y_R_F = f4.Histo1D("k_th_y_R_F")


# fill vertex histograms

# Line 1220 (k.vtx_x)
#           (k.vtx_x_L)
#           (k.vtx_x_R)
h_vtx_x    = f4.Histo1D("k_vtx_x")
h_vtx_x_L  = f4.Histo1D("k_vtx_x_L")
h_vtx_x_R  = f4.Histo1D("k_vtx_x_R")

# Line 1224 (k.vtx_y)
#           (k.vtx_y_L)
#           (k.vtx_y_R)
h_vtx_y    = f4.Histo1D("k_vtx_y")
h_vtx_y_L  = f4.Histo1D("k_vtx_y_L")
h_vtx_y_R  = f4.Histo1D("k_vtx_y_R")

# Line 1228:
#            (k.vtx_x_R, k.vtx_x_L)
#            (k.vtx_y_R, k.vtx_y_L)
h_vtx_x_L_vs_vtx_x_R = f4.Histo1D("k_vtx_x_R", "k_vtx_x_L")
h_vtx_y_L_vs_vtx_y_R = f4.Histo1D("k_vtx_y_R", "k_vtx_y_L")

# Line 1231:
#            (k.th_x_L, k.vtx_x_L)
#            (k.th_x_R, k.vtx_x_R)
#            (k.th_y_L, k.vtx_y_L)
#            (k.th_y_R, k.vtx_y_R)
h_vtx_x_L_vs_th_x_L = f4.Histo1D("k_th_x_L", "k_vtx_x_L")
h_vtx_x_R_vs_th_x_R = f4.Histo1D("k_th_x_R", "k_vtx_x_R")
h_vtx_y_L_vs_th_y_L = f4.Histo1D("k_th_y_L", "k_vtx_y_L")
h_vtx_y_R_vs_th_y_R = f4.Histo1D("k_th_y_R", "k_vtx_y_R")

# Line 1236:
#           (k.vtx_x_R - k.vtx_x_L)
#           (k.vtx_y_R - k.vtx_y_L)
h_vtx_x_diffLR = f4.Histo1D("k_vtx_x_diffLR");
h_vtx_y_diffLR = f4.Histo1D("k_vtx_y_diffLR");

# Line 1239:
#           (k.th_x, k.vtx_x_R - k.vtx_x_L)
#           (k.th_y, k.vtx_y_R - k.vtx_y_L)
h_vtx_x_diffLR_vs_th_x = f4.Histo1D("k_th_x", "k_vtx_x_diffLR");
h_vtx_y_diffLR_vs_th_y = f4.Histo1D("k_th_y", "k_vtx_y_diffLR");

# TODO TProfile
#p_vtx_x_diffLR_vs_th_x = f4.Profile1D("k_th_x", "k_vtx_y_diffLR");
#p_vtx_y_diffLR_vs_th_y = f4.Profile1D("k_th_y", "k_vtx_y_diffLR");

# Line 1245:
#           (k.vtx_x_R, k.vtx_x_R - k.vtx_x_L)
#           (k.vtx_y_R, k.vtx_y_R - k.vtx_y_L)
h_vtx_x_diffLR_vs_vtx_x_R = f4.Histo1D("k_vtx_x_R", "k_vtx_y_diffLR");
h_vtx_y_diffLR_vs_vtx_y_R = f4.Histo1D("k_vtx_y_R", "k_vtx_y_diffLR");

# TODO: safe
# if (safe)
# {
# 	h_vtx_x_safe->Fill(k.vtx_x);
# 	h_vtx_y_safe->Fill(k.vtx_y);
#
# 	h_vtx_x_diffLR_safe->Fill(k.vtx_x_R - k.vtx_x_L);
# 	h_vtx_y_diffLR_safe->Fill(k.vtx_y_R - k.vtx_y_L);
#
# 	h_vtx_x_diffLR_safe_corr->Fill((k.vtx_x_R - k.vtx_x_L) - 1080.*k.th_x);
# 	h_vtx_y_diffLR_safe_corr->Fill((k.vtx_y_R - k.vtx_y_L) + 156. *k.th_y);
# }

# TODO: min
# Line 1338
#th_y_L_min = min(th_y_L_min, k.th_y_L);
#th_y_R_min = min(th_y_R_min, k.th_y_R);


# TODO: from line 1358 to 1404 (end of event loop)

# Line 1401
# calculate acceptance divergence correction
r7 = f4.Define("correction", "CalculateAcceptanceCorrectionsRDF( kinematics )") \
       .Define("corr", "correction.corr") \
       .Define("div_corr", "correction.div_corr")

# TODO Line 1406
# for (unsigned int bi = 0; bi < binnings.size(); bi++)
# {
# 	bh_t_Nev_before[bi]->Fill(k.t, 1.);
# 	bh_t_before[bi]->Fill(k.t, 1.);
# }

# Line 1412
modelreal = ROOT.TH2D("h_th_y_vs_th_x_before", ";#theta_{x};#theta_{y}", 150, -300E-6, +300E-6, 150, -150E-6, +150E-6)
modelreal.Sumw2()
model = ROOT.RDF.TH2DModel(modelreal)
# FIXME: Weight value is missing (1.)
h_th_y_vs_th_x_before = r7.Histo2D(model, "k_th_x", "k_th_y");

# Line 1414
# Filter skip
f5 = r7.Filter("! correction.skip", "acceptance correction")

# Line 1429
modelreal = ROOT.TH1D("h_t_after", ";|t|",128, 0., 0.)
modelreal.Sumw2()
model = ROOT.RDF.TH1DModel(modelreal)
bh_t_after_ob_1_30_02 = f5.Histo1D(model, "k_t", "corr");

# Line 1435
modelreal = ROOT.TH2D("h_th_y_vs_th_x_after", ";#theta_{x};#theta_{y}", 150, -300E-6, +300E-6, 150, -150E-6, +150E-6);
modelreal.Sumw2()
model = ROOT.RDF.TH2DModel(modelreal)
h_th_y_vs_th_x_after = f5.Histo2D(model, "k_th_x", "k_th_y", "div_corr");

# Line 1435
modelreal = ROOT.TH2D("h_th_vs_phi_after", ";#phi;#theta", 50, -M_PI, +M_PI, 50, 150E-6, 550E-6);
modelreal.Sumw2()
model = ROOT.RDF.TH2DModel(modelreal)
h_th_vs_phi_after = f5.Histo2D(model, "k_th_x", "k_th_y", "div_corr");

# Line 1441
# apply normalization
modelreal = ROOT.TH1D("h_t_normalized", ";|t|",128, 0., 0.)
bh_t_normalized_ob_1_30_02 = f5.Define("corr_norm", "corr * normalization") \
                               .Histo1D("k_t", "corr_norm")

# Line 1445
modelreal = ROOT.TH2D("h_th_y_vs_th_x_normalized", ";#theta_{x};#theta_{y}", 150, -600E-6, +600E-6, 150, -600E-6, +600E-6);
modelreal.Sumw2()
model = ROOT.RDF.TH2DModel(modelreal)
h_th_y_vs_th_x_normalized = f5.Define("div_corr_norm", "correction.div_corr * normalization") \
                              .Histo2D(model, "k_th_x", "k_th_y", "div_corr_norm");


###############################
###### END OF EVENT LOOP ######
###############################


print("---------------------------- after event loop ---------------------------\n");

print(">> th_min = %s\n" % th_min);
print(">> th_y_L_min = %s\n" % th_y_L_min);
print(">> th_y_R_min = %s\n" % th_y_R_min);

print("\n");
print("N_anal = %s\n" % N_anal);
print("N_anal_zeroBias = %s\n" % N_anal_zeroBias);
print("N_zeroBias_el = %s\n" % N_zeroBias_el);
print("N_zeroBias_el_RP_trig = %s\n" % N_zeroBias_el_RP_trig);

print("N_4outof4 = %s\n" % N_4outof4);
print("N_el = %s\n" % N_el);
print("N_el_T2trig = %s\n" % N_el_T2trig);
print("N_4outof4_T2trig = %s\n" % N_4outof4_T2trig);


# TODO: From line 1454 on, mostly pure root

# Line 1492
# Normalize histograms
bh_t_after_ob_1_30_02.Scale(1., "width");

# Line 1494
bh_t_normalized_ob_1_30_02.Scale(1., "width")

# TODO Line 1497
h_th_y_vs_th_x_normalized.Scale(1., "width")

# Line 1506
th_y_diffLR.Scale(1., "width");
th_x_diffLR.Scale(1., "width");

# TODO
#th_y_diffLR_safe.Scale(1., "width");
#th_x_diffLR_safe.Scale(1., "width");

# TODO
# for (unsigned int bi = 0; bi < binnings.size(); bi++)
# 	HideLowTBins(bh_t_normalized[bi], anal.t_min_fit);

# TODO
# fit histograms
# From 1515 - 1607 (TProfiles)

# TODO
# th_y_diffLR_safe->Fit("gaus");
# th_x_diffLR_safe->Fit("gaus");

# save histograms
c = ROOT.TCanvas();

outF = ROOT.TFile.Open(outputDir+"/distributions_" + selected_diagonal + ".root", "recreate");
ROOT.gDirectory = outF.mkdir("metadata");

c = ROOT.TCanvas("rate cmp");
h_timestamp_dgn.Draw();
#h_timestamp_B0.SetLineColor(4);
#h_timestamp_B0.Draw("sames");
h_timestamp_sel.SetLineColor(2);
h_timestamp_sel.Draw("sames");
c.Write();

# Line 1700
hitDistDir = outF.mkdir("hit distributions");
ROOT.gDirectory = hitDistDir.mkdir("vertical, aligned, before selection");
h_y_L_2_F_vs_x_L_2_F_al_nosel.Write();
h_y_L_2_N_vs_x_L_2_N_al_nosel.Write();
h_y_L_1_F_vs_x_L_1_F_al_nosel.Write();
h_y_R_1_F_vs_x_R_1_F_al_nosel.Write();
h_y_R_2_N_vs_x_R_2_N_al_nosel.Write();
h_y_R_2_F_vs_x_R_2_F_al_nosel.Write();

ROOT.gDirectory = hitDistDir.mkdir("vertical, not aligned, after selection");
h_y_L_2_F_vs_x_L_2_F_noal_sel.Write();
h_y_L_2_N_vs_x_L_2_N_noal_sel.Write();
h_y_L_1_F_vs_x_L_1_F_noal_sel.Write();
h_y_R_1_F_vs_x_R_1_F_noal_sel.Write();
h_y_R_2_N_vs_x_R_2_N_noal_sel.Write();
h_y_R_2_F_vs_x_R_2_F_noal_sel.Write();

ROOT.gDirectory = hitDistDir.mkdir("vertical, aligned, after selection");
h_y_L_2_F_vs_x_L_2_F_al_sel.Write();
h_y_L_2_N_vs_x_L_2_N_al_sel.Write();
h_y_L_1_F_vs_x_L_1_F_al_sel.Write();
h_y_R_1_F_vs_x_R_1_F_al_sel.Write();
h_y_R_2_N_vs_x_R_2_N_al_sel.Write();
h_y_R_2_F_vs_x_R_2_F_al_sel.Write();

ROOT.gDirectory = outF.mkdir("selected - hits");
ROOT.gDirectory = outF.mkdir("selected - angles");

th_x_diffLR.Write()
th_y_diffLR.Write()

th_x_diffLF.Write()
th_x_diffRF.Write()

h_th_x_diffLR_vs_th_x.Write();
h_th_y_diffLR_vs_th_y.Write();
h_th_x_diffLR_vs_vtx_x.Write();

# TODO TProfile
# p_th_x_diffLR_vs_th_x->Write();
# p_th_y_diffLR_vs_th_y->Write();
# p_th_y_L_diffNF_vs_th_y_L->Write();
# p_th_y_R_diffNF_vs_th_y_R->Write();

# p_th_x_diffLR_vs_vtx_x->Write();

## TODO:
#th_x_sigmaLR_vs_th_x.Write();
#th_y_sigmaLR_vs_th_y.Write();

#th_x_diffLR_safe.Write();
#th_y_diffLR_safe.Write();

# TODO TProfile
# p_th_x_vs_th_y->Write();
# p_th_x_L_vs_th_y_L->Write();
# p_th_x_R_vs_th_y_R->Write();

h_th_y_L_vs_th_x_L.Write();
h_th_y_R_vs_th_x_R.Write();
h_th_y_vs_th_x.Write();

# TODO TProfile
# g_th_y_L_vs_th_x_L->Write();
# g_th_y_R_vs_th_x_R->Write();
# g_th_y_vs_th_x->Write();

c = ROOT.TCanvas();
c.SetLogz(1);
c.ToggleEventStatus();
c.SetCrosshair(1);
h_th_y_L_vs_th_y_R.Draw("colz");
#g_th_y_L_vs_th_y_R.Draw("p");
c.Write("canvas_th_y_L_vs_th_y_R");

h_th_x.Write();
h_th_y.Write();
h_th_y_flipped.Write();

h_th_x_L.Write();
h_th_x_R.Write();

h_th_y_L.Write();
h_th_y_R.Write();

h_th_y_L_F.Write();
h_th_y_L_N.Write();
h_th_y_R_N.Write();
h_th_y_R_F.Write();

# TODO TGraph
# {
# 	double x[] = {0, 1, 2, 3};
# 	double y[] = {anal.th_y_lcut_L, anal.th_y_hcut_L, anal.th_y_lcut_R, anal.th_y_hcut_R};
# 	TGraph *g = new TGraph(4, x, y);
# 	g->SetName("g_th_y_cuts");
# 	g->Write();
# }

ROOT.gDirectory = outF.mkdir("selected - vertex");
h_vtx_x.Write();
h_vtx_x_L.Write();
h_vtx_x_R.Write();

h_vtx_y.Write();
h_vtx_y_L.Write();
h_vtx_y_R.Write();

# # TODO:
#h_vtx_x_safe.Write();
#h_vtx_y_safe.Write();

h_vtx_x_L_vs_vtx_x_R.Write();
h_vtx_y_L_vs_vtx_y_R.Write();

h_vtx_x_L_vs_th_x_L.Write();
h_vtx_x_R_vs_th_x_R.Write();
h_vtx_y_L_vs_th_y_L.Write();
h_vtx_y_R_vs_th_y_R.Write();

h_vtx_x_diffLR.Write();
h_vtx_y_diffLR.Write();

# TODO
#h_vtx_x_diffLR_safe.Write();
#h_vtx_y_diffLR_safe.Write();

# TODO
#h_vtx_x_diffLR_safe_corr.Write();
#h_vtx_y_diffLR_safe_corr.Write();

h_vtx_x_diffLR_vs_th_x.Write();
h_vtx_y_diffLR_vs_th_y.Write();

# TODO TProfile
# p_vtx_x_diffLR_vs_th_x.Write();
# p_vtx_y_diffLR_vs_th_y.Write();

h_vtx_x_diffLR_vs_vtx_x_R.Write();
h_vtx_y_diffLR_vs_vtx_y_R.Write();


ROOT.gDirectory = outF.mkdir("optics");
#h_x_L_F_vs_th_x_L.Write();
#h_x_R_F_vs_th_x_R.Write();

# TODO
#p_x_L_N_vs_th_x_L.Write();
#p_x_L_F_vs_th_x_L.Write();
#p_x_R_N_vs_th_x_R.Write();
#p_x_R_F_vs_th_x_R.Write();

# TODO TProfile
# p_th_x_R_vs_th_x_L.Write();
# p_th_y_R_vs_th_y_L.Write();
# p_th_y_LF_vs_th_y_LN.Write();
# p_th_y_RF_vs_th_y_RN.Write();

#h_thl_y_L_vs_y_LF.Write();
#h_thl_y_R_vs_y_RF.Write();
#h_thl_y_L_vs_thl_y_R.Write();

# TODO TProfile
# p_thl_y_L_vs_y_LF.Write();
# p_thl_y_R_vs_y_RF.Write();
# p_thl_y_L_vs_thl_y_R.Write();

# TODO
# ROOT.gDirectory = outF->mkdir("binning");
# for (unsigned int bi = 0; bi < binnings.size(); bi++)
# {
# 	TGraph *g = new TGraph();
# 	g->SetName(("g_binning_"+binnings[bi]).c_str());
# 	g->SetTitle(";bin center;bin width");
#
# 	TH1D *h = bh_t_Nev_before[bi];
# 	for (int bin = 1; bin <= h->GetNbinsX(); bin++)
# 		g->SetPoint(g->GetN(), h->GetBinCenter(bin), h->GetBinWidth(bin));
#
# 	g->Write();
# }

# TODO TGraph/TProfile
# ROOT.gDirectory = outF->mkdir("time dependences");
# p_diffLR_th_x_vs_time->Write();
# ProfileToRMSGraph(p_diffLR_th_x_vs_time, gRMS_diffLR_th_x_vs_time);
# gRMS_diffLR_th_x_vs_time->Write();
#
# p_diffLR_th_y_vs_time->Write();
# ProfileToRMSGraph(p_diffLR_th_y_vs_time, gRMS_diffLR_th_y_vs_time);
# gRMS_diffLR_th_y_vs_time->Write();
#
# p_diffNF_th_y_L_vs_time->Write();
# ProfileToRMSGraph(p_diffNF_th_y_L_vs_time, gRMS_diffNF_th_y_L_vs_time);
# gRMS_diffNF_th_y_L_vs_time->Write();
#
# p_diffNF_th_y_R_vs_time->Write();
# ProfileToRMSGraph(p_diffNF_th_y_R_vs_time, gRMS_diffNF_th_y_R_vs_time);
# gRMS_diffNF_th_y_R_vs_time->Write();
#
# p_vtx_x_vs_time->Write();
# ProfileToRMSGraph(p_vtx_x_vs_time, gRMS_vtx_x_vs_time);
# gRMS_vtx_x_vs_time->Write();

# TGraphErrors *g_beam_div_x_vs_time = new TGraphErrors; g_beam_div_x_vs_time->SetName("g_beam_div_x_vs_time"); g_beam_div_x_vs_time->SetTitle(";timestamp;beam divergence in x");
# TGraphErrors *g_sensor_res_x_vs_time = new TGraphErrors; g_sensor_res_x_vs_time->SetName("g_sensor_res_x_vs_time"); g_sensor_res_x_vs_time->SetTitle(";timestamp;sensor resolution in x");
# for (int i = 0; i <= gRMS_vtx_x_vs_time->GetN(); ++i)
# {
# 	double time=0., si_diff=0., si_vtx=0.;
# 	gRMS_vtx_x_vs_time->GetPoint(i, time, si_vtx);
# 	gRMS_diffLR_th_x_vs_time->GetPoint(i, time, si_diff);
#
# 	double si_bdx = si_vtx * sqrt(2.) / 90. * 1E-3;	// in rad
# 	double si_srx = sqrt(si_diff*si_diff/2. - si_bdx*si_bdx);
#
# 	g_beam_div_x_vs_time->SetPoint(i, time, si_bdx);
# 	g_sensor_res_x_vs_time->SetPoint(i, time, si_srx);
# }
#
# g_beam_div_x_vs_time->Write();
# g_sensor_res_x_vs_time->Write();
#
# p_th_x_R_vs_time->Write();
# p_th_y_R_vs_time->Write();
# p_th_x_L_vs_time->Write();
# p_th_y_L_vs_time->Write();
#
# g_ext_diffLR_th_x_vs_time->Write();
# g_ext_diffLR_th_y_vs_time->Write();
#
# p_input_beam_div_x_vs_time->Write();
# p_input_beam_div_y_vs_time->Write();
#
# g_L_L_F_vs_time->Write();
# g_L_R_F_vs_time->Write();
#

accDir = outF.mkdir("acceptance correction");
ROOT.gDirectory = accDir.mkdir("ob_1_30_02");
bh_t_after_ob_1_30_02.Write()
# for (unsigned int bi = 0; bi < binnings.size(); bi++)
# {
# 	ROOT.gDirectory = accDir->mkdir(binnings[bi].c_str());
# 	bh_t_Nev_before[bi]->Write();
# 	bh_t_Nev_after_no_corr[bi]->Write();
# 	bh_t_before[bi]->Write();
# 	bh_t_after_no_corr[bi]->Write();
# 	bh_t_after[bi]->Write();
# 	bp_t_phi_corr[bi]->Write();
# 	bp_t_full_corr[bi]->Write();
#
# 	c = new TCanvas("t cmp");
# 	c->SetLogy(1);
# 	bh_t_after[bi]->Draw("");
# 	bh_t_before[bi]->Draw("same");
# 	c->Write();
# }

# FIXME accDir should consider binnigs
ROOT.gDirectory = accDir;

# TODO
# p_t_ub_div_corr->Write();

# h_th_y_vs_th_x_before.Write();
# h_th_y_vs_th_x_after.Write();
# h_th_vs_phi_after.Write();

# g_weight_vs_th_y->Write();
#
# g_th_y_vs_th_x_acc->Write();
#
normDir = outF.mkdir("normalization");
ROOT.gDirectory = normDir.mkdir("ob_1_30_02");
bh_t_normalized_ob_1_30_02.Write()
# for (unsigned int bi = 0; bi < binnings.size(); bi++)
# {
# 	ROOT.gDirectory = normDir->mkdir(binnings[bi].c_str());
# 	bh_t_normalized[bi]->Write();
# 	bh_t_normalized_rel_diff[bi]->Write();
# }
#
# p_norm_corr->Write();
# p_3outof4_corr->Write();
#
# h_th_y_vs_th_x_normalized->Write();
#
# g_norm_corr_vs_div_corr->Write();
#
# TDirectory *normUnfDir = outF->mkdir("normalization+unfolding");
# for (unsigned int bi = 0; bi < binnings.size(); bi++)
# {
# 	ROOT.gDirectory = normUnfDir->mkdir(binnings[bi].c_str());
#
# 	bh_t_normalized_unsmeared[bi]->Write();
# 	bh_t_normalized_unsmeared_rel_diff[bi]->Write();
# }
