import ROOT
import sys
import os.path

RDF = ROOT.ROOT.RDataFrame

# Load C++ headers
ROOT.gInterpreter.Declare('#include "common_definitions.h"')
ROOT.gInterpreter.Declare('#include "parameters_global.h"')
ROOT.gInterpreter.Declare('#include "common_algorithms.h"')

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

# Read all branches
rdf = RDF(treename, input_file)

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
initAnalysis="""
// select cuts
extern Analysis anal;
anal.BuildCuts();
anal.n_si = input_n_si;
"""
ROOT.gInterpreter.Declare(initAnalysis)

# Line 241 - Let's start assuming no overrideCutSelection
# if (overrideCutSelection)

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

# TODO Line 380
# get unsmearing correction

# Line 394
# book metadata histograms
# timestamp_bins = timestamp_max - timestamp_min + 1.;

# Long TODO
# Lines 397 - 779
# THistograms, TGraphs and TProfiles declaration

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
# // remove troublesome runs
# unsigned int run = ev.run_num / 100000;
# unsigned int file = ev.run_num % 100000;
# if (SkipRun(run, file, true))
# 	continue;
#
# // update timestamp run boundaries
# auto rtbit = runTimestampBoundaries.find(run);
# if (rtbit == runTimestampBoundaries.end())
# {
# 	runTimestampBoundaries.insert({run, {ev.timestamp, ev.timestamp}});
# } else {
# 	rtbit->second.first = min(rtbit->second.first, ev.timestamp);
# 	rtbit->second.second = max(rtbit->second.second, ev.timestamp);
# }
#
# // check time - selected?
# if (anal.SkipTime(ev.timestamp))
# 	continue;
#
# if (time_group_divisor != 0)
# {
# 	double time_group_interval = 1.;	// s
# 	int time_group = int(ev.timestamp / time_group_interval);
# 	if ( (time_group % time_group_divisor) != time_group_remainder)
# 		continue;
# }

# Diagonal cut (L831)
f1 = rdf.Filter("v_L_2_F && v_L_2_N && v_R_2_F && v_R_2_N", 'allDiagonalRPs')

# FIXME: N_4outof4_T2trig is not used to produce any plot
# Line 835
# if ((ev.trigger_bits & 64) != 0)
# 	N_4outof4_T2trig++;

# Line 838
h_timestamp_dgn = f1.Histo1D("timestamp")

# Line 841
# select the elastic-trigger bunch(es) only
# if (SkipBunch(run, ev.bunch_num))
#
# SkipBunch will return False for whatever DS in use,
# since all of them define keepAllBunches = True within
# parameters.h

# zero bias event?
# Create named filter with number of zero bias events
f2 = f1.Filter("(trigger_bits & 512) != 0", 'zero_bias_event')

# Apply fine alignment (L 852)


# run reconstruction (Line 868)
xs = ["x_L_1_F", "x_L_2_N", "x_L_2_F", "x_R_1_F", "x_R_2_N", "x_R_2_F"]
ys = ["y_L_1_F", "y_L_2_N", "y_L_2_F", "y_R_1_F", "y_R_2_N", "y_R_2_F"]

## apply fine alignment
r2 = f1.Define("h_al", "ApplyFineAlignment( timestamp ,{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {})".format(*(xs+ys) )) \
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
# h_al = ev.h
# so filling a hist with h_al.L_1_F.x equals to use the branch x_L_1_F
h_y_L_1_F_vs_x_L_1_F_al_nosel = r2.Histo1D("h_al_x_L_1_F", "h_al_y_L_1_F")
h_y_L_2_N_vs_x_L_2_N_al_nosel = r2.Histo1D("h_al_x_L_2_N", "h_al_y_L_2_N")
h_y_L_2_F_vs_x_L_2_F_al_nosel = r2.Histo1D("h_al_x_L_2_F", "h_al_y_L_2_F")
h_y_R_1_F_vs_x_R_1_F_al_nosel = r2.Histo1D("h_al_x_R_1_F", "h_al_y_R_1_F")
h_y_R_2_N_vs_x_R_2_N_al_nosel = r2.Histo1D("h_al_x_R_2_N", "h_al_y_R_2_N")
h_y_R_2_F_vs_x_R_2_F_al_nosel = r2.Histo1D("h_al_x_R_2_F", "h_al_y_R_2_F")

# TODO
#if (detailsLevel >= 2)
#{
#    h_timestamp_B0->Fill(ev.timestamp);
#    g_run_vs_timestamp->SetPoint(g_run_vs_timestamp->GetN(), ev.timestamp, ev.run_num);
#    //g_ev_num_vs_timestamp->SetPoint(g_ev_num_vs_timestamp->GetN(), ev.timestamp, ev.event_num);
#    //g_tr_num_vs_timestamp->SetPoint(g_tr_num_vs_timestamp->GetN(), ev.timestamp, ev.trigger_num);
#}


## Line 877
### kinematics struct
ks = r2.Define("kinematics", 'DoReconstruction( h_al )')

ks_ext= ks.Define("k_th_x_R",        "kinematics.th_x_R") \
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

## TODO fill no-cut histograms (L957)
# for (unsigned int ci = 1; ci <= anal.N_cuts; ++ci)
# {
# 	//h2_cq_full[ci]->Fill(ccb[ci]*cqa[ci] - cca[ci]*cqb[ci], cca[ci]*cqa[ci] + ccb[ci]*cqb[ci] + ccc[ci]);
# 	h2_cq_full[ci]->Fill(cd.cqa[ci], cd.cqb[ci]);
# }

# TODO (L988)
# g_selected_bunch_num_vs_timestamp->SetPoint(g_selected_bunch_num_vs_timestamp->GetN(), ev.timestamp, ev.bunch_num);

# TODO (L1040)
# p_norm_corr->Fill(ev.timestamp, norm_corr);
# p_3outof4_corr->Fill(k.th_y, inefficiency_3outof4);


# TODO:
# // data for alignment
# // (SHOULD use hit positions WITHOUT alignment corrections, i.e. ev.h)
# signed int period = int((ev.timestamp - anal.alignment_t0) / anal.alignment_ts);

# TODO fill histograms at line 1110
h_y_L_1_F_vs_x_L_1_F_noal_sel = f1.Histo1D("x_L_1_F", "y_L_1_F")
h_y_L_2_N_vs_x_L_2_N_noal_sel = f1.Histo1D("x_L_2_N", "y_L_2_N")
h_y_L_2_F_vs_x_L_2_F_noal_sel = f1.Histo1D("x_L_2_F", "y_L_2_F")
h_y_R_1_F_vs_x_R_1_F_noal_sel = f1.Histo1D("x_R_1_F", "y_R_1_F")
h_y_R_2_N_vs_x_R_2_N_noal_sel = f1.Histo1D("x_R_2_N", "y_R_2_N")
h_y_R_2_F_vs_x_R_2_F_noal_sel = f1.Histo1D("x_R_2_F", "y_R_2_F")

h_y_L_1_F_vs_x_L_1_F_al_sel = f1.Histo1D("x_L_1_F", "y_L_1_F")
h_y_L_2_N_vs_x_L_2_N_al_sel = f1.Histo1D("x_L_2_N", "y_L_2_N")
h_y_L_2_F_vs_x_L_2_F_al_sel = f1.Histo1D("x_L_2_F", "y_L_2_F")
h_y_R_1_F_vs_x_R_1_F_al_sel = f1.Histo1D("x_R_1_F", "y_R_1_F")
h_y_R_2_N_vs_x_R_2_N_al_sel = f1.Histo1D("x_R_2_N", "y_R_2_N")
h_y_R_2_F_vs_x_R_2_F_al_sel = f1.Histo1D("x_R_2_F", "y_R_2_F")

# Line 1157 (k.th_x_R - k.th_x_L)
#           (k.th_y_R - k.th_y_L)
th_x_diffLR = ks_ext.Histo1D("k_th_x_diffLR")
th_y_diffLR = ks_ext.Histo1D("k_th_y_diffLR")

# Line 1160 (k.th_x_L - k.th_x)
#           (k.th_x_R - k.th_x)
th_x_diffLF = ks_ext.Histo1D("k_th_x_diffLF")
th_x_diffRF = ks_ext.Histo1D("k_th_x_diffRF")

# Line 1163 (k.th_x, k.th_x_R - k.th_x_L)
#           (k.th_y, k.th_y_R - k.th_y_L)
#           (k.vtx_x, k.th_x_R - k.th_x_L)
h_th_x_diffLR_vs_th_x  = ks_ext.Histo1D("k_th_x", "k_th_x_diffLR")
h_th_y_diffLR_vs_th_y  = ks_ext.Histo1D("k_th_y", "k_th_y_diffLR")
h_th_x_diffLR_vs_vtx_x = ks_ext.Histo1D("k_vtx_x", "k_th_x_diffLR")

# TODO:
# Line 1168 TProfile
#p_th_x_diffLR_vs_th_x = ks_ext.Histo1D("k_th_x", "k_th_x_diffLR")
#p_th_y_diffLR_vs_th_y = ks_ext.Histo1D("k_th_y", "k_th_y_diffLR")
#p_th_y_L_diffNF_vs_th_y_L = ks_ext.Histo1D("k_th_y_L", "k_th_y_L_diffNF")
#p_th_y_R_diffNF_vs_th_y_R = ks_ext.Histo1D("k_th_y_R", "k_th_y_R_diffNF")

# TODO:
# Line 1172 TProfile
#p_th_x_diffLR_vs_vtx_x = ks_ext.Histo1D("k_vtx_x", "k_th_x_diffRF")

# TODO:
# if (safe)
# {
# 	th_y_diffLR_safe->Fill(k.th_y_R - k.th_y_L);
# 	th_x_diffLR_safe->Fill(k.th_x_R - k.th_x_L);
# }

# TODO:
# Line 1184 TProfile
#p_th_x_vs_th_y = ks_ext.Histo1D("k_th_y", "k_th_x")
#p_th_x_L_vs_th_y_L = ks_ext.Histo1D("k_th_y_L", "k_th_x_L")
#p_th_x_R_vs_th_y_R = ks_ext.Histo1D("k_th_y_R", "k_th_x_R")

# Line 1188 (k.th_x_L, k.th_y_L)
#           (k.th_x_R, k.th_y_R)
#           (k.th_x, k.th_y)
h_th_y_L_vs_th_x_L = ks_ext.Histo1D("k_th_x_L", "k_th_y_L")
h_th_y_R_vs_th_x_R = ks_ext.Histo1D("k_th_x_R", "k_th_y_R")
h_th_y_vs_th_x     = ks_ext.Histo1D("k_th_x", "k_th_y")

# TODO: TGraph
# Line 1192
# if (detailsLevel >= 1)
# {
# 	g_th_y_L_vs_th_x_L->SetPoint(g_th_y_L_vs_th_x_L->GetN(), k.th_x_L, k.th_y_L);
# 	g_th_y_R_vs_th_x_R->SetPoint(g_th_y_R_vs_th_x_R->GetN(), k.th_x_R, k.th_y_R);
# 	g_th_y_vs_th_x->SetPoint(g_th_y_vs_th_x->GetN(), k.th_x, k.th_y);
# }

# Line 1199 (k.th_y_R, k.th_y_L)
h_th_y_L_vs_th_y_R = ks_ext.Histo1D("k_th_y_R", "k_th_y_L")

# TODO TGraph
# if (detailsLevel > 2)
# 	g_th_y_L_vs_th_y_R->SetPoint(g_th_y_L_vs_th_y_R->GetN(), k.th_y_R, k.th_y_L);

# Line 1203: (k.th_x)
#            (k.th_y)
h_th_x     = ks_ext.Histo1D("k_th_x")
h_th_y     = ks_ext.Histo1D("k_th_y")

# Line 1205: (-k.th_y)
h_th_y_flipped = ks_ext.Histo1D("minus_k_th_y")

# Line 1207: (k.th_x_L)
#            (k.th_x_R)
h_th_x_L   = ks_ext.Histo1D("k_th_x_L")
h_th_x_R   = ks_ext.Histo1D("k_th_x_R")

# Line 1210: (k.th_y_L)
#            (k.th_y_R)
h_th_y_L   = ks_ext.Histo1D("k_th_y_L")
h_th_y_R   = ks_ext.Histo1D("k_th_y_R")

# Line 1213: (k.th_y_L_F)
#            (k.th_y_L_N)
#            (k.th_y_R_N)
#            (k.th_y_R_F)
h_th_y_L_F = ks_ext.Histo1D("k_th_y_L_F")
h_th_y_L_N = ks_ext.Histo1D("k_th_y_L_N")
h_th_y_R_N = ks_ext.Histo1D("k_th_y_R_N")
h_th_y_R_F = ks_ext.Histo1D("k_th_y_R_F")


# fill vertex histograms

# Line 1220 (k.vtx_x)
#           (k.vtx_x_L)
#           (k.vtx_x_R)
h_vtx_x    = ks_ext.Histo1D("k_vtx_x")
h_vtx_x_L  = ks_ext.Histo1D("k_vtx_x_L")
h_vtx_x_R  = ks_ext.Histo1D("k_vtx_x_R")

# Line 1224 (k.vtx_y)
#           (k.vtx_y_L)
#           (k.vtx_y_R)
h_vtx_y    = ks_ext.Histo1D("k_vtx_y")
h_vtx_y_L  = ks_ext.Histo1D("k_vtx_y_L")
h_vtx_y_R  = ks_ext.Histo1D("k_vtx_y_R")

# Line 1228:
#            (k.vtx_x_R, k.vtx_x_L)
#            (k.vtx_y_R, k.vtx_y_L)
h_vtx_x_L_vs_vtx_x_R = ks_ext.Histo1D("k_vtx_x_R", "k_vtx_x_L")
h_vtx_y_L_vs_vtx_y_R = ks_ext.Histo1D("k_vtx_y_R", "k_vtx_y_L")

# Line 1231:
#            (k.th_x_L, k.vtx_x_L)
#            (k.th_x_R, k.vtx_x_R)
#            (k.th_y_L, k.vtx_y_L)
#            (k.th_y_R, k.vtx_y_R)
h_vtx_x_L_vs_th_x_L = ks_ext.Histo1D("k_th_x_L", "k_vtx_x_L")
h_vtx_x_R_vs_th_x_R = ks_ext.Histo1D("k_th_x_R", "k_vtx_x_R")
h_vtx_y_L_vs_th_y_L = ks_ext.Histo1D("k_th_y_L", "k_vtx_y_L")
h_vtx_y_R_vs_th_y_R = ks_ext.Histo1D("k_th_y_R", "k_vtx_y_R")

# Line 1236:
#           (k.vtx_x_R - k.vtx_x_L)
#           (k.vtx_y_R - k.vtx_y_L)
h_vtx_x_diffLR = ks_ext.Histo1D("k_vtx_x_diffLR");
h_vtx_y_diffLR = ks_ext.Histo1D("k_vtx_y_diffLR");

# Line 1239:
#           (k.th_x, k.vtx_x_R - k.vtx_x_L)
#           (k.th_y, k.vtx_y_R - k.vtx_y_L)
h_vtx_x_diffLR_vs_th_x = ks_ext.Histo1D("k_th_x", "k_vtx_x_diffLR");
h_vtx_y_diffLR_vs_th_y = ks_ext.Histo1D("k_th_y", "k_vtx_y_diffLR");

# TODO TProfile
#p_vtx_x_diffLR_vs_th_x = ks_ext.Profile1D("k_th_x", "k_vtx_y_diffLR");
#p_vtx_y_diffLR_vs_th_y = ks_ext.Profile1D("k_th_y", "k_vtx_y_diffLR");

# Line 1245:
#           (k.vtx_x_R, k.vtx_x_R - k.vtx_x_L)
#           (k.vtx_y_R, k.vtx_y_R - k.vtx_y_L)
h_vtx_x_diffLR_vs_vtx_x_R = ks_ext.Histo1D("k_vtx_x_R", "k_vtx_y_diffLR");
h_vtx_y_diffLR_vs_vtx_y_R = ks_ext.Histo1D("k_vtx_y_R", "k_vtx_y_diffLR");

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


# TODO: from line 1358 to 1452 (end of event loop)

# TODO: From line 1454 on, mostly pure root

# TODO Line 1497
# h_th_y_vs_th_x_normalized.Scale(1., "width")

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
# TODO timestamp
#h_timestamp_B0.SetLineColor(4);
#h_timestamp_B0.Draw("sames");
#h_timestamp_sel.SetLineColor(2);
#h_timestamp_sel.Draw("sames");
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
# TDirectory *normDir = outF->mkdir("normalization");
# for (unsigned int bi = 0; bi < binnings.size(); bi++)
# {
# 	ROOT.gDirectory = normDir->mkdir(binnings[bi].c_str());
# 	bh_t_normalized[bi]->Write();
# 	bh_t_normalized_rel_diff[bi]->Write();
# }
#
# ROOT.gDirectory = normDir;
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
