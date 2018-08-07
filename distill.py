import ROOT
import glob
import sys

RDF = ROOT.ROOT.RDataFrame

# Branches clasified by diagonal
diagonals = {
    # return a tuple: ([left] verticals in 45, [right] verticals in 56))
    "d45b_56t"  : (["track_rp_5", "track_rp_21", "track_rp_25"],
                   ["track_rp_104", "track_rp_120", "track_rp_124"]),
    "ad45b_56b" : (["track_rp_5", "track_rp_21", "track_rp_25"],
                   ["track_rp_105", "track_rp_121", "track_rp_125"]),
    "d45t_56b"  : (["track_rp_4", "track_rp_20", "track_rp_24"],
                   ["track_rp_105", "track_rp_121", "track_rp_125"]),
    "ad45t_56t" : (["track_rp_4", "track_rp_20", "track_rp_24"],
                   ["track_rp_104", "track_rp_120", "track_rp_124"])
}

DS = {
    'DS1'   : '4495', 
    'DS2'   : '4496', 
    'DS3'   : '4499', 
    'DS4'   : '4505', 
    'DS5'   : '4509', 
    'DS6'   : '4510', 
    'DS7'   : '4511'
}

threads_description = "no_MT"

if len(sys.argv) < 3:
    print('Usage: python distill.py <diagonal> <DS> [threads number]')
    sys.exit(1)  # no diagonal specified

if len(sys.argv) == 4:
    if int(sys.argv[3]) < 1:
        print('Threads number should be > 0')
        sys.exit(1)  # wrong threads number
    ROOT.ROOT.EnableImplicitMT(int(sys.argv[3]))
    threads_description = "threads_" + sys.argv[3]

# Select branches
selected_diagonal = sys.argv[1]
selected_DS = sys.argv[2]
if selected_diagonal not in diagonals.keys():
    print('Invalid diagonal: %s' % selected_diagonal)
    print('Choose between:   %s' % diagonals.keys())
    sys.exit(1)

if selected_DS not in DS.keys():
    print('DS not available: %s' % selected_DS)
    print('Choose between: %s' % DS.keys())
    sys.exit(1)

rp_left, rp_right = diagonals[selected_diagonal]

# Extracted from: DS1/block1/input_files.h
source_file       = "input_files_{}.txt".format(selected_DS)
input_ntuple_name = "TotemNtuple"
prefix            = "root://eostotem.cern.ch//eos/totem/data/cmstotem/2015/90m/Totem/Ntuple/version2/{}/".format(DS[selected_DS])
input_files       = [prefix + line.rstrip('\n') for line in open(source_file)]

# Convert to PyROOT vector
vec_input_files = ROOT.vector('string')()
[vec_input_files.push_back(f) for f in input_files]

# Columns per branch
attributes = ['valid', 'x', 'y']

full_branches = ["{}.{}".format(c,a) for a in attributes for c in rp_left+rp_right ]

# Split left and right branch on valid, x and y
valids = [ "(unsigned int) {}".format(v) for v in full_branches[0:6]]
xs     = full_branches[6:12]
ys     = full_branches[12:18]

print("Selected branches: \n" + "\n\t".join(full_branches))

# Filter and define output branches
filter_code = """({0}.valid + {1}.valid + {2}.valid ) >= 2 &&
({3}.valid + {4}.valid + {5}.valid ) >= 2
""".format(*(rp_left+rp_right))

print("Filter code: \n" + filter_code)

# Input tree
treename= "TotemNtuple"
rdf = RDF(treename, vec_input_files)

# Output tree, file and branches
outTreeName = "distilled"
outFileName = "distill_{}_{}_{}_new.root".format(selected_DS, threads_description, selected_diagonal)
branchList  = ["v_L_1_F", "x_L_1_F", "y_L_1_F",
               "v_L_2_N", "x_L_2_N", "y_L_2_N",
               "v_L_2_F", "x_L_2_F", "y_L_2_F",
               "v_R_1_F", "x_R_1_F", "y_R_1_F",
               "v_R_2_N", "x_R_2_N", "y_R_2_N",
               "v_R_2_F", "x_R_2_F", "y_R_2_F",
               "timestamp",
               "run_num",
               "bunch_num",
               "event_num",
               "trigger_num",
               "trigger_bits"
              ]

# Convert to PyROOT vector
vec_outbranchlist = ROOT.vector('string')()
[vec_outbranchlist.push_back(b) for b in branchList]

# Filter and define output branches
r = rdf.Filter(filter_code)  \
       .Define("v_L_1_F", valids[0]) \
       .Define("x_L_1_F", xs[0]) \
       .Define("y_L_1_F", ys[0]) \
       .Define("v_L_2_N", valids[1]) \
       .Define("x_L_2_N", xs[1]) \
       .Define("y_L_2_N", ys[1]) \
       .Define("v_L_2_F", valids[2]) \
       .Define("x_L_2_F", xs[2]) \
       .Define("y_L_2_F", ys[2]) \
       .Define("v_R_1_F", valids[3]) \
       .Define("x_R_1_F", xs[3]) \
       .Define("y_R_1_F", ys[3]) \
       .Define("v_R_2_N", valids[4]) \
       .Define("x_R_2_N", xs[4]) \
       .Define("y_R_2_N", ys[4]) \
       .Define("v_R_2_F", valids[5]) \
       .Define("x_R_2_F", xs[5]) \
       .Define("y_R_2_F", ys[5]) \
       .Define("timestamp",    "(unsigned int) (event_info.timestamp - 1444860000)") \
       .Define("run_num",      "(unsigned int) event_info.run_no")                 \
       .Define("bunch_num",    "trigger_data.bunch_num")            \
       .Define("event_num",    "trigger_data.event_num")            \
       .Define("trigger_num",  "trigger_data.trigger_num")          \
       .Define("trigger_bits", "trigger_data.input_status_bits")

# All above actions are not executed at the moment they are called,
# but they are lazy, i.e. delayed until the moment one of their results
# is accessed (in this case by .GetValue() )
print("Distilled events: %s" % r.Count().GetValue())

# Save output tree
r.Snapshot(outTreeName, outFileName, vec_outbranchlist)
