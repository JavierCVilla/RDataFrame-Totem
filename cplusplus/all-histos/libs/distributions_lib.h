#ifndef _distributions_lib_h_
#define _distributions_lib_h_

#include <string>
#include <vector>
#include <set>
#include <map>
#include <cmath>

#include "TMath.h"
#include "TGraph.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDSymEigen.h"
#include "TRandom2.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

//**********************
// common_definitions.h
//**********************


enum DiagonalType { dUnknown, d45b_56t, d45t_56b, dCombined, ad45b_56b, ad45t_56t };

DiagonalType diagonal = dUnknown;

double th_y_sign = 0.;

//----------------------------------------------------------------------------------------------------

struct AlignmentData
{
    //    a: xy coupling in rad
    //    b: x shift in mm
    //    c: y shift in mm
    double a_L_2_F, b_L_2_F, c_L_2_F;
    double a_L_2_N, b_L_2_N, c_L_2_N;
    double a_L_1_F, b_L_1_F, c_L_1_F;

    double a_R_1_F, b_R_1_F, c_R_1_F;
    double a_R_2_N, b_R_2_N, c_R_2_N;
    double a_R_2_F, b_R_2_F, c_R_2_F;

    AlignmentData();

    /*
    AlignmentData Interpolate(double s_N, double s_F, double s_NH, double s_FH) const
    {
        AlignmentData r;

        r.a_L_F = a_L_N + (a_L_F - a_L_N)/(s_F - s_N) * (s_FH - s_N); r.a_L_N = a_L_N + (a_L_F - a_L_N)/(s_F - s_N) * (s_NH - s_N);
        r.a_R_F = a_R_N + (a_R_F - a_R_N)/(s_F - s_N) * (s_FH - s_N); r.a_R_N = a_R_N + (a_R_F - a_R_N)/(s_F - s_N) * (s_NH - s_N);

        r.b_L_F = b_L_N + (b_L_F - b_L_N)/(s_F - s_N) * (s_FH - s_N); r.b_L_N = b_L_N + (b_L_F - b_L_N)/(s_F - s_N) * (s_NH - s_N);
        r.b_R_F = b_R_N + (b_R_F - b_R_N)/(s_F - s_N) * (s_FH - s_N); r.b_R_N = b_R_N + (b_R_F - b_R_N)/(s_F - s_N) * (s_NH - s_N);

        r.c_L_F = c_L_N + (c_L_F - c_L_N)/(s_F - s_N) * (s_FH - s_N); r.c_L_N = c_L_N + (c_L_F - c_L_N)/(s_F - s_N) * (s_NH - s_N);
        r.c_R_F = c_R_N + (c_R_F - c_R_N)/(s_F - s_N) * (s_FH - s_N); r.c_R_N = c_R_N + (c_R_F - c_R_N)/(s_F - s_N) * (s_NH - s_N);

        return r;
    }
    */
};

//----------------------------------------------------------------------------------------------------

enum AlignmentType { atNone, atConstant, atTimeDependent };

struct AlignmentSource
{
    struct GraphSet
    {
        TGraph *L_2_F, *L_2_N, *L_1_F, *R_1_F, *R_2_N, *R_2_F;
        GraphSet() : L_2_F(NULL), L_2_N(NULL), L_1_F(NULL), R_1_F(NULL), R_2_N(NULL), R_2_F(NULL)
        {
        }
    } gs_a, gs_b, gs_c;

    AlignmentData cnst;

    AlignmentType type_a, type_b, type_c;
    string src_a, src_b, src_c;

    AlignmentSource() : type_a(atNone), type_b(atNone), type_c(atNone)
    {
    }

    void SetAlignmentA(AlignmentType t, const string &fn = "");
    void SetAlignmentB(AlignmentType t, const string &fn = "");
    void SetAlignmentC(AlignmentType t, const string &fn = "");

    void InitOne(const string label, AlignmentType t, const string &fn, GraphSet &gs, const string &obj);

    void Init();

    AlignmentData Eval(double timestamp) const;
};

//----------------------------------------------------------------------------------------------------

struct UnitHitData
{
    // validity flag
    unsigned int v;

    // hit position in mm
    double x, y;

    UnitHitData() : v(0), x(0.), y(0.) {}

    void operator += (const UnitHitData &add);
};

//----------------------------------------------------------------------------------------------------

struct HitData
{
    UnitHitData L_1_F, L_2_N, L_2_F;
    UnitHitData R_1_F, R_2_N, R_2_F;


    void operator += (const HitData &add);

    HitData ApplyAlignment(const AlignmentData &al) const;

    // TODO: remove hard-coded z positions
    /*
    HitData ApplyInterpolatedAlignment(const AlignmentData &a, double sN, double sF) const
    {
        AlignmentData a_int = a.Interpolate(214.628, 220.000, sN, sF);

        return ApplyAlignment(a_int);
    }
    */
};

//----------------------------------------------------------------------------------------------------

struct EventRed
{
    unsigned int timestamp;
    unsigned int run_num, bunch_num, event_num, trigger_num;
    unsigned int trigger_bits;

    // vertical RPs
    HitData h;

    //HitData hH;    // horizontal RPs
};

//----------------------------------------------------------------------------------------------------

struct Environment
{
    // beam momentum (GeV)
    double p, p_L, p_R;

    // beam momentum uncertainty
    double si_de_p;

    // beam divergence
    double si_th_x_L, si_th_y_L;        // rad
    double si_th_x_R, si_th_y_R;        // rad

    double si_th_y_RL_assym_unc;        // uncertainty of the L-R assymetry

    // vertex smearing
    double si_vtx_x, si_vtx_y;         // mm

    // pitch-induced error
    double si_de_P_L, si_de_P_R;    // mm

    // optics
    double v_x_L_1_F, v_x_L_2_N, v_x_L_2_F, v_x_R_1_F, v_x_R_2_N, v_x_R_2_F;    // 1
    double v_y_L_1_F, v_y_L_2_N, v_y_L_2_F, v_y_R_1_F, v_y_R_2_N, v_y_R_2_F;    // 1
    double L_x_L_1_F, L_x_L_2_N, L_x_L_2_F, L_x_R_1_F, L_x_R_2_N, L_x_R_2_F;    // mm
    double L_y_L_1_F, L_y_L_2_N, L_y_L_2_F, L_y_R_1_F, L_y_R_2_N, L_y_R_2_F;    // mm

    // optics: x-y coupling (x = L_x * th_x + v_x * x^* + la_x * th_y)
    /*
    double la_x_L_F, la_x_L_N, la_x_R_N, la_x_R_F;    // mm
    double la_y_L_F, la_y_L_N, la_y_R_N, la_y_R_F;    // mm
    */

    // optics perturbation covariance matrices
    // order of elements:
    //        left arm:  v_x_L_N, L_x_L_N, v_y_L_N, L_y_L_N, v_x_L_F, L_x_L_F, v_y_L_F, L_y_L_F
    //        right arm: v_x_R_N, L_x_R_N, v_y_R_N, L_y_R_N, v_x_R_F, L_x_R_F, v_y_R_F, L_y_R_F
    // units: v's in 1, L's in m
    TMatrixDSym opt_cov;

    // optics perturbation generator matrices
    TMatrixD opt_per_gen;

    // alignment uncertainties
    double si_de_x, si_de_y_R, si_de_y_D, si_tilt;

    // misalignments (mm)
    double de_x_L_N, de_y_L_N, tilt_L_N;
    double de_x_L_F, de_y_L_F, tilt_L_F;
    double de_x_R_N, de_y_R_N, tilt_R_N;
    double de_x_R_F, de_y_R_F, tilt_R_F;

    Environment() : opt_cov(16), opt_per_gen(16, 16)
    {
    }

    void InitNominal();
    void UseMatchedOptics();

    void PrintOpticsUncertainties() const;

    void Print() const;

    void ApplyRandomOpticsPerturbations(TVectorD &de);

    void ApplyRandomOpticsPerturbations();

    /// modes counted from 0 to 15
    void ApplyOpticsPerturbationMode(int mode, double coef);

    /// modes counted from 0 to 7
    void ApplyEffectiveLengthPerturbationMode(int mode, double coef);
};

//----------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------

struct Kinematics
{
    double th_x_L_F, th_x_L_N, th_x_R_N, th_x_R_F, th_x_L, th_x_R, th_x;    //    rad
    double th_y_L_F, th_y_L_N, th_y_R_N, th_y_R_F, th_y_L, th_y_R, th_y;    //    rad

    double vtx_x_L_F, vtx_x_L_N, vtx_x_R_N, vtx_x_R_F, vtx_x_L, vtx_x_R, vtx_x;    // in mm
    double vtx_y_L_F, vtx_y_L_N, vtx_y_R_N, vtx_y_R_F, vtx_y_L, vtx_y_R, vtx_y;    // in mm

    double th;                // in rad
    double phi;                // in rad
    double t_x, t_y, t;        // in GeV^2

    Kinematics() : th_y(0.) {}

    void ThetasToTPhi(const Environment &env);

    void TPhiToThetas(const Environment &env);
};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

struct CutData
{
    double cqa[9];    ///< array of quantities qa
    double cqb[9];    ///< array of quantities qb
    double cv[9];    ///< array of cut quantities v = a*qa + b*qb + c
    bool ct[9];        ///< array of flags whether |v| < n_si * si
    bool select;

    CutData() :  select(true) {}
};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

struct Analysis
{
    // binning, |t| in GeV^2
    double t_min, t_max;
    double t_min_full, t_max_full;
    double t_min_fit;

    // elastic selection cuts
    double n_si;

    double cut1_a, cut1_c, cut1_si;
    double cut2_a, cut2_c, cut2_si;
    double cut3_a, cut3_c, cut3_si;
    double cut4_a, cut4_c, cut4_si;
    double cut5_a, cut5_c, cut5_si;
    double cut6_a, cut6_c, cut6_si;
    double cut7_a, cut7_c, cut7_si;
    double cut8_a, cut8_c, cut8_si;

    std::vector< std::pair<double, double> > timeIntervals;

    unsigned int N_cuts;    // number of cuts - indexed from 1!
    string cqaN[9], cqbN[9];
    double cca[9], ccb[9], ccc[9], csi[9];
    std::vector<unsigned int> cuts;    // list of active cuts

    // analysis cuts (rad)
    double th_y_lcut_L, th_y_lcut_R, th_y_lcut;
    double th_y_hcut_L, th_y_hcut_R, th_y_hcut;

    double th_x_lcut;
    double th_x_hcut;

    // (un)-smearing parameters
    double si_th_x_1arm_L;
    double si_th_x_1arm_R;
    double si_th_x_1arm_unc;
    double si_th_x_2arm;
    double si_th_x_2arm_unc;

    double si_th_y_1arm;
    double si_th_y_1arm_unc;
    double si_th_y_2arm;
    double si_th_y_2arm_unc;

    // efficiency parameters
    bool use_3outof4_efficiency_fits;        // whether to use time-dependent fits of 3-out-of-4 efficiency
    bool use_pileup_efficiency_fits;        // whether to use time-dependent fits of pile-up efficiency

    double inefficiency_3outof4;            // inefficiency from 3-out-of-4 method, used only if use_3outof4_efficiency_fits=false
    double inefficiency_shower_near;        // inefficiency due to shower in near RP
    double inefficiency_pile_up;            // inefficiency due to pile-up, used only if use_pileup_efficiency_fits=false
    double inefficiency_trigger;            // trigger inefficiency
    double inefficiency_DAQ;                // DAQ inefficiency

    // normalisation correction to subtract background
    double bckg_corr;

    // (delivered) luminosity
    double L_int;    // mb^-1

    // 3-out-of-4 efficiency uncertainty (only used in MC simulation)
    double eff_3outof4_fixed_point, eff_3outof4_slope, eff_3outof4_slope_unc;

    // normalisation correction and its uncertainty (only used in MC simulation)
    double norm_corr, norm_corr_unc;

    double alignment_t0;    // beginning of the first time-slice
    double alignment_ts;    // time-slice in s

    double eff_th_y_min;

    // y ranges for alignment
    struct AlignmentYRange
    {
        double bot_min, bot_max, top_min, top_max;
        AlignmentYRange(double bmi=0., double bma=0., double tmi=0., double tma=0.) :
            bot_min(bmi), bot_max(bma), top_min(tmi), top_max(tma) {}
    };
    map<std::string, AlignmentYRange> alignmentYRanges;

    void BuildCuts();
    bool EvaluateCuts(const HitData &, const Kinematics &, CutData &) const;

    bool SkipTime(unsigned int timestamp) const;

    void Print() const;
};

CutData EvaluateCutsRDF( const HitData &h_al, const Kinematics &k );

//----------------------------------------------------------------------------------------------------

struct Correction
{
    double phi_corr;
    double div_corr;
    double corr;
    bool skip;

    Correction() : phi_corr(0.), div_corr(0.), corr(0.), skip(true) {}
};

//----------------------------------------------------------------------------------------------------

struct Binning {
    unsigned int N_bins;
    double *bin_edges;
};

//********************** End common_definitions ******************

//----------------------------------------------------------------------------------------------------

//**********************
// parameters_global.h
//**********************

double timestamp0 = 1444860000;

string storageDir;

vector<string> distilledNtuples;

vector<AlignmentSource> alignmentSources;
Analysis anal;
Environment env;

string unsmearing_file;
string unsmearing_object;

string luminosity_data_file;

//----------------------------------------------------------------------------------------------------

void Init_global();

//----------------------------------------------------------------------------------------------------

void Init_global_45b_56t();

//----------------------------------------------------------------------------------------------------

void Init_global_45t_56b();

//----------------------------------------------------------------------------------------------------

//********************** End parameters_global ******************

//----------------------------------------------------------------------------------------------------

//**********************
// common_algorithms.h
//**********************

//----------------------------------------------------------------------------------------------------

// Old Kinematics struct has been modified (end of this file)

//----------------------------------------------------------------------------------------------------

void BuildBinning(const Analysis &anal, const string &type, double* &binEdges, unsigned int &bins, bool verbose = false);

void BuildBinningRDF(const Analysis &anal, const string &type, Binning &b);

//----------------------------------------------------------------------------------------------------

bool CalculateAcceptanceCorrections(double th_y_sign,
        const Kinematics &k, const Analysis &anal,
        double &phi_corr, double &div_corr);

Correction CalculateAcceptanceCorrectionsRDF( const Kinematics &k);

//----------------------------------------------------------------------------------------------------

bool SkipRun(unsigned int /*run*/, unsigned int /*file*/, bool /*strict = true */);

//----------------------------------------------------------------------------------------------------

// map: run number (8372) --> list of triggered bunches
typedef std::map<unsigned int, std::vector<unsigned int> > BunchMap;

bool keepAllBunches;
BunchMap bunchMap;

bool SkipBunch(unsigned int run, unsigned bunch);

//----------------------------------------------------------------------------------------------------

// returns the beam for which the bunch is non-colliding
// for colliding bunches returns zero
unsigned int NonCollidingBunch(unsigned int /*run*/, unsigned /*bunch*/);

//----------------------------------------------------------------------------------------------------

bool IsZeroBias(unsigned int trigger, unsigned int /*run*/, unsigned int /*event*/);

//----------------------------------------------------------------------------------------------------

HitData ProtonTransport(const Kinematics & /*k*/, const Environment & /*env*/);

HitData ApplyFineAlignment( unsigned int &timestamp,
                            double &x_L_1_F, double &x_L_2_N, double &x_L_2_F,
                            double &x_R_1_F, double &x_R_2_N, double &x_R_2_F,
                            double &y_L_1_F, double &y_L_2_N, double &y_L_2_F,
                            double &y_R_1_F, double &y_R_2_N, double &y_R_2_F);

Kinematics DoReconstruction(HitData &h);

// Wrapper around anal.Skiptime
bool SkipTime( unsigned int &timestamp);

// Custom function to replace original check in line distributions.cc::820
bool SkipTimeInterval( unsigned int &timestamp, int &tgd, int &tgr );

// Custom function to replace original check in line distributions.cc::1021
double getNorm_corr( unsigned int &timestamp );

// Custom function to replace original check in line distributions.cc::1048
double getNormalization( double &norm_corr );

// FIXME Optimize this
// This functions is meant to be used in a RDF::Define
// where a column will be defined containing a 1 value for event
double One();

//********************** End common_algorithms ******************

//----------------------------------------------------------------------------------------------------

//**********************
// parameters.h
//**********************

double timestamp_min = 20.9E3, timestamp_max = 31.5E3;

void Init_base();

//----------------------------------------------------------------------------------------------------

void Init_45b_56t();

//----------------------------------------------------------------------------------------------------

void Init_45t_56b();

//********************** End parameters_global ******************

//----------------------------------------------------------------------------------------------------

//**********************
// common.h
//**********************

int rcIncompatibleDiagonal = 123;

//----------------------------------------------------------------------------------------------------

void Init(const std::string &dgnStr);

#endif

