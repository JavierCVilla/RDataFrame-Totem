#include "distributions_lib.h"

/******************************************************************
 * common_definitions.h
 ******************************************************************/

using namespace std;

AlignmentData::AlignmentData()
{
    a_L_2_F = b_L_2_F = c_L_2_F = 0.;
    a_L_2_N = b_L_2_N = c_L_2_N = 0.;
    a_L_1_F = b_L_1_F = c_L_1_F = 0.;

    a_R_1_F = b_R_1_F = c_R_1_F = 0.;
    a_R_2_N = b_R_2_N = c_R_2_N = 0.;
    a_R_2_F = b_R_2_F = c_R_2_F = 0.;
}

void AlignmentSource::SetAlignmentA(AlignmentType t, const string &fn)
{
    type_a = t;
    src_a = fn;
}

void AlignmentSource::SetAlignmentB(AlignmentType t, const string &fn)
{
    type_b = t;
    src_b = fn;
}

void AlignmentSource::SetAlignmentC(AlignmentType t, const string &fn)
{
    type_c = t;
    src_c = fn;
}

void AlignmentSource::InitOne(const string label, AlignmentType t, const string &fn, GraphSet &gs, const string &obj)
{
    printf(">> AlignmentSource::InitOne > alignment `%s': type %u\n", label.c_str(), t);

    if (t == atTimeDependent)
    {
        TFile *alF = new TFile(fn.c_str());

        if (alF->IsZombie())
        {
            printf("\tERROR: cannot open file with alignment graphs.\n");
            delete alF;
            return;
        }

        TGraph *g_L_2_F = (TGraph *) alF->Get(( string("L_2_F/") + obj).c_str() );
        TGraph *g_L_2_N = (TGraph *) alF->Get(( string("L_2_N/") + obj).c_str() );
        TGraph *g_L_1_F = (TGraph *) alF->Get(( string("L_1_F/") + obj).c_str() );

        TGraph *g_R_1_F = (TGraph *) alF->Get(( string("R_1_F/") + obj).c_str() );
        TGraph *g_R_2_N = (TGraph *) alF->Get(( string("R_2_N/") + obj).c_str() );
        TGraph *g_R_2_F = (TGraph *) alF->Get(( string("R_2_F/") + obj).c_str() );

        if (g_L_2_F && g_L_2_N && g_L_1_F && g_R_1_F && g_R_2_N && g_R_2_F)
        {
            printf("\talignment graphs successfully loaded\n");
        } else {
            printf("\tERROR: unable to load some alignment graphs\n");
            delete alF;
            return;
        }

        gs.L_2_F = new TGraph(*g_L_2_F);
        gs.L_2_N = new TGraph(*g_L_2_N);
        gs.L_1_F = new TGraph(*g_L_1_F);

        gs.R_1_F = new TGraph(*g_R_1_F);
        gs.R_2_N = new TGraph(*g_R_2_N);
        gs.R_2_F = new TGraph(*g_R_2_F);

        delete alF;
    }
}

void AlignmentSource::Init()
{
    printf(">> AlignmentSource::Init\n");
    InitOne("a", type_a, src_a, gs_a, "a_fit");
    InitOne("b", type_b, src_b, gs_b, "b_fit");
    InitOne("c", type_c, src_c, gs_c, "c_fit");
}

AlignmentData AlignmentSource::Eval(double timestamp) const
{
    AlignmentData d;

    // a
    if (type_a == atNone)
    {
        d.a_L_2_F = 0.;
        d.a_L_2_N = 0.;
        d.a_L_1_F = 0.;

        d.a_R_1_F = 0.;
        d.a_R_2_N = 0.;
        d.a_R_2_F = 0.;
    }

    if (type_a == atConstant)
    {
        d.a_L_2_F = cnst.a_L_2_F;
        d.a_L_2_N = cnst.a_L_2_N;
        d.a_L_1_F = cnst.a_L_1_F;

        d.a_R_1_F = cnst.a_R_1_F;
        d.a_R_2_N = cnst.a_R_2_N;
        d.a_R_2_F = cnst.a_R_2_F;
    }

    if (type_a == atTimeDependent)
    {
        d.a_L_2_F = gs_a.L_2_F->Eval(timestamp)*1E-3;
        d.a_L_2_N = gs_a.L_2_N->Eval(timestamp)*1E-3;
        d.a_L_1_F = gs_a.L_1_F->Eval(timestamp)*1E-3;

        d.a_R_1_F = gs_a.R_1_F->Eval(timestamp)*1E-3;
        d.a_R_2_N = gs_a.R_2_N->Eval(timestamp)*1E-3;
        d.a_R_2_F = gs_a.R_2_F->Eval(timestamp)*1E-3;
    }

    // b
    if (type_b == atNone)
    {
        d.b_L_2_F = 0.;
        d.b_L_2_N = 0.;
        d.b_L_1_F = 0.;

        d.b_R_1_F = 0.;
        d.b_R_2_N = 0.;
        d.b_R_2_F = 0.;
    }

    if (type_b == atConstant)
    {
        d.b_L_2_F = cnst.b_L_2_F;
        d.b_L_2_N = cnst.b_L_2_N;
        d.b_L_1_F = cnst.b_L_1_F;

        d.b_R_1_F = cnst.b_R_1_F;
        d.b_R_2_N = cnst.b_R_2_N;
        d.b_R_2_F = cnst.b_R_2_F;
    }

    if (type_b == atTimeDependent)
    {
        d.b_L_2_F = gs_b.L_2_F->Eval(timestamp)*1E-3;
        d.b_L_2_N = gs_b.L_2_N->Eval(timestamp)*1E-3;
        d.b_L_1_F = gs_b.L_1_F->Eval(timestamp)*1E-3;

        d.b_R_1_F = gs_b.R_1_F->Eval(timestamp)*1E-3;
        d.b_R_2_N = gs_b.R_2_N->Eval(timestamp)*1E-3;
        d.b_R_2_F = gs_b.R_2_F->Eval(timestamp)*1E-3;
    }

    // c
    if (type_c == atNone)
    {
        d.c_L_2_F = 0.;
        d.c_L_2_N = 0.;
        d.c_L_1_F = 0.;

        d.c_R_1_F = 0.;
        d.c_R_2_N = 0.;
        d.c_R_2_F = 0.;
    }

    if (type_c == atConstant)
    {
        d.c_L_2_F = cnst.c_L_2_F;
        d.c_L_2_N = cnst.c_L_2_N;
        d.c_L_1_F = cnst.c_L_1_F;

        d.c_R_1_F = cnst.c_R_1_F;
        d.c_R_2_N = cnst.c_R_2_N;
        d.c_R_2_F = cnst.c_R_2_F;
    }

    if (type_c == atTimeDependent)
    {
        d.c_L_2_F = gs_c.L_2_F->Eval(timestamp)*1E-3;
        d.c_L_2_N = gs_c.L_2_N->Eval(timestamp)*1E-3;
        d.c_L_1_F = gs_c.L_1_F->Eval(timestamp)*1E-3;

        d.c_R_1_F = gs_c.R_1_F->Eval(timestamp)*1E-3;
        d.c_R_2_N = gs_c.R_2_N->Eval(timestamp)*1E-3;
        d.c_R_2_F = gs_c.R_2_F->Eval(timestamp)*1E-3;
    }

    return d;
}

void UnitHitData::operator += (const UnitHitData &add)
{
    x += add.x;
    y += add.y;
}

void HitData::operator += (const HitData &add)
{
    L_1_F += add.L_1_F;
    L_2_N += add.L_2_N;
    L_2_F += add.L_2_F;

    R_1_F += add.R_1_F;
    R_2_N += add.R_2_N;
    R_2_F += add.R_2_F;
}

HitData HitData::ApplyAlignment(const AlignmentData &al) const
{
    HitData r;

    r.L_2_F.x = L_2_F.x - al.a_L_2_F * L_2_F.y - al.b_L_2_F; r.L_2_F.y = L_2_F.y - al.c_L_2_F;
    r.L_2_N.x = L_2_N.x - al.a_L_2_N * L_2_N.y - al.b_L_2_N; r.L_2_N.y = L_2_N.y - al.c_L_2_N;
    r.L_1_F.x = L_1_F.x - al.a_L_1_F * L_1_F.y - al.b_L_1_F; r.L_1_F.y = L_1_F.y - al.c_L_1_F;

    r.R_1_F.x = R_1_F.x - al.a_R_1_F * R_1_F.y - al.b_R_1_F; r.R_1_F.y = R_1_F.y - al.c_R_1_F;
    r.R_2_N.x = R_2_N.x - al.a_R_2_N * R_2_N.y - al.b_R_2_N; r.R_2_N.y = R_2_N.y - al.c_R_2_N;
    r.R_2_F.x = R_2_F.x - al.a_R_2_F * R_2_F.y - al.b_R_2_F; r.R_2_F.y = R_2_F.y - al.c_R_2_F;

    return r;
}

void Environment::Print() const
{
    printf("p=%E, p_L=%E, p_R=%E\n", p, p_L, p_R);
    printf("\n");
    printf("si_th_x_L=%E, si_th_y_L=%E\n", si_th_x_L, si_th_y_L);
    printf("si_th_x_R=%E, si_th_y_R=%E\n", si_th_x_R, si_th_y_R);
    printf("si_vtx_x=%E, si_vtx_y=%E\n", si_vtx_x, si_vtx_y);
    printf("si_de_P_L=%E, si_de_P_R=%E\n", si_de_P_L, si_de_P_R);
    printf("\n");

    printf("L_x_L_1_F = %E, v_x_L_1_F = %E, L_y_L_1_F = %E, v_y_L_1_F = %E\n", L_x_L_1_F, v_x_L_1_F, L_y_L_1_F, v_y_L_1_F);
    printf("L_x_L_2_N = %E, v_x_L_2_N = %E, L_y_L_2_N = %E, v_y_L_2_N = %E\n", L_x_L_2_N, v_x_L_2_N, L_y_L_2_N, v_y_L_2_N);
    printf("L_x_L_2_F = %E, v_x_L_2_F = %E, L_y_L_2_F = %E, v_y_L_2_F = %E\n", L_x_L_2_F, v_x_L_2_F, L_y_L_2_F, v_y_L_2_F);
    printf("L_x_R_1_F = %E, v_x_R_1_F = %E, L_y_R_1_F = %E, v_y_R_1_F = %E\n", L_x_R_1_F, v_x_R_1_F, L_y_R_1_F, v_y_R_1_F);
    printf("L_x_R_2_N = %E, v_x_R_2_N = %E, L_y_R_2_N = %E, v_y_R_2_N = %E\n", L_x_R_2_N, v_x_R_2_N, L_y_R_2_N, v_y_R_2_N);
    printf("L_x_R_2_F = %E, v_x_R_2_F = %E, L_y_R_2_F = %E, v_y_R_2_F = %E\n", L_x_R_2_F, v_x_R_2_F, L_y_R_2_F, v_y_R_2_F);

    printf("\n");
    printf("si_de_x=%E, si_de_y_R=%E, si_de_y_D=%E, si_tilt=%E\n", si_de_x, si_de_y_R, si_de_y_D, si_tilt);
    printf("\n");
    printf("de_x_L_N=%E, de_y_L_N=%E, tilt_L_N=%E\n", de_x_L_N, de_y_L_N, tilt_L_N);
    printf("de_x_L_F=%E, de_y_L_F=%E, tilt_L_F=%E\n", de_x_L_F, de_y_L_F, tilt_L_F);
    printf("de_x_R_N=%E, de_y_R_N=%E, tilt_R_N=%E\n", de_x_R_N, de_y_R_N, tilt_R_N);
    printf("de_x_R_F=%E, de_y_R_F=%E, tilt_R_F=%E\n", de_x_R_F, de_y_R_F, tilt_R_F);
    printf("\n");
    printf("si_th_y_RL_assym_unc=%E\n", si_th_y_RL_assym_unc);

    PrintOpticsUncertainties();
}

void Environment::ApplyRandomOpticsPerturbations()
{
    TVectorD de(16);
    ApplyRandomOpticsPerturbations(de);
}

void Environment::ApplyRandomOpticsPerturbations(TVectorD & /*de*/){};


//----------------------------------------------------------------------------------------------------

void Environment::ApplyOpticsPerturbationMode(int /*mode*/, double /*coef*/)
{
    /*
    printf(">> Environment::ApplyOpticsPerturbationMode\n");

    // prepare correlation matrix
    TMatrixDSym cor(opt_cov);
    TMatrixDSym Sigma(opt_cov);
    for (int i = 0; i < opt_cov.GetNrows(); i++)
        for (int j = 0; j < opt_cov.GetNcols(); j++)
        {
            cor(i, j) /= sqrt( opt_cov(i, i) * opt_cov(j, j) );
            Sigma(i, j) = (i == j) ? sqrt( opt_cov(i, i) ) : 0.;
        }

    // eigen decomposition
    TMatrixDSymEigen eig_decomp(cor);
    TVectorD eig_values(eig_decomp.GetEigenValues());

    // construct mode
    TVectorD vm(opt_cov.GetNrows());
    for (int i = 0; i < opt_cov.GetNrows(); i++)
    {
        double l = eig_values(i);
        double sl = (l > 0.) ? sqrt(l) : 0.;
        vm(i) = (i == mode) ? sl * coef : 0.;
    }

    vm = Sigma * eig_decomp.GetEigenVectors() * vm;

    printf("\tleft arm: mode %u, coefficient %+.3f\n", mode, coef);
    vm.Print();

    v_x_L_N += vm(0) * 1E0;
    L_x_L_N += vm(1) * 1E3;
    v_y_L_N += vm(2) * 1E0;
    L_y_L_N += vm(3) * 1E3;
    v_x_L_F += vm(4) * 1E0;
    L_x_L_F += vm(5) * 1E3;
    v_y_L_F += vm(6) * 1E0;
    L_y_L_F += vm(7) * 1E3;

    v_x_R_N += vm(8) * 1E0;
    L_x_R_N += vm(9) * 1E3;
    v_y_R_N += vm(10) * 1E0;
    L_y_R_N += vm(11) * 1E3;
    v_x_R_F += vm(12) * 1E0;
    L_x_R_F += vm(13) * 1E3;
    v_y_R_F += vm(14) * 1E0;
    L_y_R_F += vm(15) * 1E3;
    */
}

//----------------------------------------------------------------------------------------------------

void Environment::ApplyEffectiveLengthPerturbationMode(int /*mode*/, double /*coef*/)
{
    /*

    printf(">> Environment::ApplyEffectiveLengthPerturbationMode\n");

    // prepare reduced covariance matrix
    TMatrixDSym cov_red(8);
    for (unsigned int i = 0; i < 8; i++)
        for (unsigned int j = 0; j < 8; j++)
            cov_red(i, j) = opt_cov(2*i+1, 2*j+1);

    // eigen decomposition
    TMatrixDSymEigen eig_decomp(cov_red);
    TVectorD eig_values(eig_decomp.GetEigenValues());

    // construct mode
    TVectorD vm(cov_red.GetNrows());
    for (int i = 0; i < cov_red.GetNrows(); i++)
    {
        double l = eig_values(i);
        double sl = (l > 0.) ? sqrt(l) : 0.;
        vm(i) = (i == mode) ? sl * coef : 0.;
    }

    vm = eig_decomp.GetEigenVectors() * vm;

    printf("\tmode %u, coefficient %+.3f\n", mode, coef);
    //vm.Print();

    L_x_L_N += vm(0) * 1E3;
    L_y_L_N += vm(1) * 1E3;
    L_x_L_F += vm(2) * 1E3;
    L_y_L_F += vm(3) * 1E3;
    L_x_R_N += vm(4) * 1E3;
    L_y_R_N += vm(5) * 1E3;
    L_x_R_F += vm(6) * 1E3;
    L_y_R_F += vm(7) * 1E3;

    */
}

//----------------------------------------------------------------------------------------------------

void Environment::PrintOpticsUncertainties() const
{
    printf("optics uncertainties: left arm\n");
    printf("\tv_x_N: %.4f\n", sqrt(opt_cov(0, 0)));
    printf("\tL_x_N: %.4f m\n", sqrt(opt_cov(1, 1)));
    printf("\tv_y_N: %.4f\n", sqrt(opt_cov(2, 2)));
    printf("\tL_y_N: %.4f m\n", sqrt(opt_cov(3, 3)));
    printf("\tv_x_F: %.4f\n", sqrt(opt_cov(4, 4)));
    printf("\tL_x_F: %.4f m\n", sqrt(opt_cov(5, 5)));
    printf("\tv_y_F: %.4f\n", sqrt(opt_cov(6, 6)));
    printf("\tL_y_F: %.4f m\n", sqrt(opt_cov(7, 7)));

    printf("optics uncertainties: right arm\n");
    printf("\tv_x_N: %.4f\n", sqrt(opt_cov(8, 8)));
    printf("\tL_x_N: %.4f m\n", sqrt(opt_cov(9, 9)));
    printf("\tv_y_N: %.4f\n", sqrt(opt_cov(10, 10)));
    printf("\tL_y_N: %.4f m\n", sqrt(opt_cov(11, 11)));
    printf("\tv_x_F: %.4f\n", sqrt(opt_cov(12, 12)));
    printf("\tL_x_F: %.4f m\n", sqrt(opt_cov(13, 13)));
    printf("\tv_y_F: %.4f\n", sqrt(opt_cov(14, 14)));
    printf("\tL_y_F: %.4f m\n", sqrt(opt_cov(15, 15)));
}

//----------------------------------------------------------------------------------------------------

void Kinematics::ThetasToTPhi(const Environment &env)
{
    th = sqrt(th_x*th_x + th_y*th_y);
    t_x = env.p*env.p * th_x * th_x;
    t_y = env.p*env.p * th_y * th_y;
    t = t_x + t_y;
    phi = atan2(th_y, th_x);
}

void Kinematics::TPhiToThetas(const Environment &env)
{
    th = sqrt(t) / env.p;
    th_x_L = th_x_R = th_x = th * cos(phi);
    th_y_L = th_y_R = th_y = th * sin(phi);

    t_x = t * cos(phi) * cos(phi);
    t_y = t * sin(phi) * sin(phi);
}


bool Analysis::SkipTime(unsigned int timestamp) const
{
    if (timeIntervals.size() == 0)
        return false;

    bool selected = false;
    for (unsigned int i = 0; i < timeIntervals.size(); i++)
    {
        if (timestamp >= timeIntervals[i].first && timestamp <= timeIntervals[i].second)
        {
            selected = true;
            break;
        }
    }

    return !selected;
}


void Analysis::Print() const
{
    printf("t_min=%E, t_max=%E, t_min_full=%E, t_max_full=%E\n", t_min, t_max, t_min_full, t_max_full);
    printf("t_min_fit=%E\n", t_min_fit);

    printf("\n");
    printf("%lu time intervals:\n", timeIntervals.size());
    for (std::vector< std::pair<double, double> >::const_iterator it = timeIntervals.begin(); it != timeIntervals.end(); ++it)
        printf("\tfrom %.1f to %.1f\n", it->first, it->second);

    printf("\n");
    printf("n_si=%E\n", n_si);

    printf("\n");
    printf("cut1_a=%E, cut1_c=%E, cut1_si=%E\n", cut1_a, cut1_c, cut1_si);
    printf("cut2_a=%E, cut2_c=%E, cut2_si=%E\n", cut2_a, cut2_c, cut2_si);
    printf("cut3_a=%E, cut3_c=%E, cut3_si=%E\n", cut3_a, cut3_c, cut3_si);
    printf("cut4_a=%E, cut4_c=%E, cut4_si=%E\n", cut4_a, cut4_c, cut4_si);
    printf("cut5_a=%E, cut5_c=%E, cut5_si=%E\n", cut5_a, cut5_c, cut5_si);
    printf("cut6_a=%E, cut6_c=%E, cut6_si=%E\n", cut6_a, cut6_c, cut6_si);
    printf("cut7_a=%E, cut7_c=%E, cut7_si=%E\n", cut7_a, cut7_c, cut7_si);
    printf("cut8_a=%E, cut8_c=%E, cut8_si=%E\n", cut8_a, cut8_c, cut8_si);

    printf("\n");
    printf("cut parameters:\n");
    for (unsigned int i = 1; i <= N_cuts; i++)
    {
        printf("%u| cqaN=%s, cqbN=%s | cca=%E, ccb=%E, ccc=%E, csi=%E\n", i,
            cqaN[i].c_str(), cqbN[i].c_str(), cca[i], ccb[i], ccc[i], csi[i]);
    }

    printf("\n");
    printf("%lu enabled cuts: ", cuts.size());
    for (unsigned int i = 0; i < cuts.size(); i++)
        printf((i == 0) ? "%i" : ", %i", cuts[i]);

    printf("\n");
    printf("th_x_lcut=%E\n", th_x_lcut);
    printf("th_x_hcut=%E\n", th_x_hcut);
    printf("th_y_lcut_L=%E, th_y_lcut_R=%E, th_y_lcut=%E\n", th_y_lcut_L, th_y_lcut_R, th_y_lcut);
    printf("th_y_hcut_L=%E, th_y_hcut_R=%E, th_y_hcut=%E\n", th_y_hcut_L, th_y_hcut_R, th_y_hcut);

    printf("\n");
    printf("si_th_x_1arm_L=%E, si_th_x_1arm_R=%E, si_th_x_1arm_unc=%E\n", si_th_x_1arm_L, si_th_x_1arm_R, si_th_x_1arm_unc);
    printf("si_th_x_2arm=%E, si_th_x_2arm_unc=%E\n", si_th_x_2arm, si_th_x_2arm_unc);
    printf("si_th_y_1arm=%E, si_th_y_1arm_unc=%E\n", si_th_y_1arm, si_th_y_1arm_unc);
    printf("si_th_y_2arm=%E, si_th_y_2arm_unc=%E\n", si_th_y_2arm, si_th_y_2arm_unc);

    printf("\n");
    printf("use_3outof4_efficiency_fits = %i\n", use_3outof4_efficiency_fits);
    printf("use_pileup_efficiency_fits= %i\n", use_pileup_efficiency_fits);
    printf("inefficiency_3outof4 = %.3f\n", inefficiency_3outof4);
    printf("inefficiency_shower_near = %.3f\n", inefficiency_shower_near);
    printf("inefficiency_pile_up = %.3f\n", inefficiency_pile_up);
    printf("inefficiency_trigger = %.3f\n", inefficiency_trigger);
    printf("inefficiency_DAQ = %.3f\n", inefficiency_DAQ);
    printf("bckg_corr = %.3f\n", bckg_corr);
    printf("L_int=%E\n", L_int);
    printf("eff_3outof4_fixed_point=%E, eff_3outof4_slope=%E, eff_3outof4_slope_unc=%E\n", eff_3outof4_fixed_point, eff_3outof4_slope, eff_3outof4_slope_unc);
    printf("norm_corr=%E, norm_corr_unc=%E\n", norm_corr, norm_corr_unc);
}

//********************** End common_definitions ******************

/******************************************************************
 * parameters_global.h
 ******************************************************************/


 void Init_global()
 {
     // environment settings
     env.InitNominal();

     // binning
     // TODO
     anal.t_min = 0.02; anal.t_max = 3.5;
     anal.t_min_full = 0.; anal.t_max_full = 4.0;

     // approximate (time independent) resolutions
     // TODO
     anal.si_th_y_1arm = 3.1E-6 / sqrt(2.);
     anal.si_th_y_1arm_unc = 0.E-6 / sqrt(2.);

     anal.si_th_y_2arm = anal.si_th_y_1arm / sqrt(2.);
     anal.si_th_y_2arm_unc = 0E-6;

     anal.si_th_x_1arm_L = 0E-6;
     anal.si_th_x_1arm_R = 0E-6;
     anal.si_th_x_1arm_unc = 0E-6;

     anal.si_th_x_2arm = 0E-6;
     anal.si_th_x_2arm_unc = 0E-6;

     // analysis settings
     anal.th_x_lcut = -1.;
     anal.th_x_hcut = +1.;
 }

 //----------------------------------------------------------------------------------------------------

 void Init_global_45b_56t()
 {
     anal.th_y_lcut_L = 30E-6; anal.th_y_lcut_R = 33.5E-6; anal.th_y_lcut = 34.5E-6;
     anal.th_y_hcut_L = 102E-6; anal.th_y_hcut_R = 102E-6; anal.th_y_hcut = 100E-6;
 }

 //----------------------------------------------------------------------------------------------------

 void Init_global_45t_56b()
 {
     anal.th_y_lcut_L = 27E-6; anal.th_y_lcut_R = 27.5E-6; anal.th_y_lcut = 28.5E-6;
     anal.th_y_hcut_L = 102E-6; anal.th_y_hcut_R = 102E-6; anal.th_y_hcut = 100E-6;
 }


/******************************************************************
 * common_cuts (defined in common_definitions.h)
 ******************************************************************/

 void Analysis::BuildCuts()
 {
     N_cuts = 8;

     // cut structure:
     //    | a*qa + b*qb + c| < n_si * si

     // a: th_x_R, b: th_x_L
     cqaN[1] = "#theta_{x}^{R}"; cqbN[1] = "#theta_{x}^{L}";
     cca[1] = -cut1_a;
     ccb[1] = 1.;
     ccc[1] = cut1_c;
     csi[1] = cut1_si;
     cuts.push_back(1);

     // a: th_y_R, b: th_y_L
     cqaN[2] = "#theta_{y}^{R}"; cqbN[2] = "#theta_{y}^{L}";
     cca[2] = -cut2_a;
     ccb[2] = 1.;
     ccc[2] = cut2_c;
     csi[2] = cut2_si;
     cuts.push_back(2);

     // a: th_x_R, b: vtx_x_R
     cqaN[3] = "#theta_{x}^{R}"; cqbN[3] = "vtx_{x}^{R}";
     cca[3] = -cut3_a;
     ccb[3] = 1.;
     ccc[3] = cut3_c;
     csi[3] = cut3_si;
     //cuts.push_back(3);

     // a: th_x_L, b: vtx_x_L
     cqaN[4] = "#theta_{x}^{L}"; cqbN[4] = "vtx_{x}^{L}";
     cca[4] = -cut4_a;
     ccb[4] = 1.;
     ccc[4] = cut4_c;
     csi[4] = cut4_si;
     //cuts.push_back(4);

     // a: y_R_N, b: y_R_F - y_R_N
     cqaN[5] = "y^{R,N}"; cqbN[5] = "y^{R,F} - y^{R,N}";
     cca[5] = -cut5_a;
     ccb[5] = 1.;
     ccc[5] = cut5_c;
     csi[5] = cut5_si;
     //cuts.push_back(5);

     // a: y_L_N, b: y_L_F - y_L_N
     cqaN[6] = "y^{L,N}"; cqbN[6] = "y^{L,F} - y^{L,N}";
     cca[6] = -cut6_a;
     ccb[6] = 1.;
     ccc[6] = cut6_c;
     csi[6] = cut6_si;
     //cuts.push_back(6);

     // a: th_x, b: vtx_x_R - vtx_x_L
     cqaN[7] = "#theta_{x}"; cqbN[7] = "vtx_{x}^{R} - vtx_{x}^{L}";
     cca[7] = -cut7_a;
     ccb[7] = 1.;
     ccc[7] = cut7_c;
     csi[7] = cut7_si;
     cuts.push_back(7);

     // a: th_y, b: vtx_y_R - vtx_y_L
     cqaN[8] = "#theta_{y}"; cqbN[8] = "vtx_{y}^{R} - vtx_{y}^{L}";
     cca[8] = -cut8_a;
     ccb[8] = 1.;
     ccc[8] = cut8_c;
     csi[8] = cut8_si;
 }

 //----------------------------------------------------------------------------------------------------

 bool Analysis::EvaluateCuts(const HitData & h, const Kinematics &k, CutData &cd) const
 {
     cd.cqa[1] = k.th_x_R;    cd.cqb[1] = k.th_x_L;
     cd.cqa[2] = k.th_y_R;    cd.cqb[2] = k.th_y_L;
     cd.cqa[3] = k.th_x_R;    cd.cqb[3] = k.vtx_x_R;
     cd.cqa[4] = k.th_x_L;    cd.cqb[4] = k.vtx_x_L;
     cd.cqa[5] = h.R_2_N.y;    cd.cqb[5] = h.R_2_F.y - h.R_2_N.y;
     cd.cqa[6] = h.L_2_N.y;    cd.cqb[6] = h.L_2_F.y - h.L_2_N.y;
     cd.cqa[7] = k.th_x;        cd.cqb[7] = k.vtx_x_R - k.vtx_x_L;
     cd.cqa[8] = k.th_y;        cd.cqb[8] = k.vtx_y_R - k.vtx_y_L;

     for (unsigned int ci = 1; ci <= N_cuts; ++ci)
     {
         cd.cv[ci] = cca[ci]*cd.cqa[ci] + ccb[ci]*cd.cqb[ci] + ccc[ci];
         cd.ct[ci] = (fabs(cd.cv[ci]) <= n_si * csi[ci]);
         //printf("cut %u: |%+E| < %E * %E <==> %i\n", ci, cd.cv[ci], n_si, csi[ci], cd.ct[ci]);
     }

     // and between all cuts
     bool select = true;
     for (unsigned int ci = 0; ci < cuts.size(); ci++)
     {
         select &= cd.ct[cuts[ci]];
     }

     return select;
 }


 //----------------------------------------------------------------------------------------------------

 CutData EvaluateCutsRDF( const HitData &h_al, const Kinematics &k ){
     CutData cd;
     extern Analysis anal;
     bool select = anal.EvaluateCuts( h_al, k, cd);
     cd.select = select;
     return cd;
 }
 /******************************************************************
 * common_parameters (no header), defined in common_definitions.h
 ******************************************************************/

 void Environment::InitNominal()
 {
     // beam momentum (GeV)
     p = p_L = p_R = 6500.;

     // momentum uncertainty
     // TODO: update
     si_de_p = 1E-3 * p;

     // angular (one-side) beam smearing (rad)
     // TODO: update
     si_th_x_L = si_th_x_R = 0.88E-6;
     si_th_y_L = si_th_y_R = 0.95E-6 / sqrt(2.);

     // vertex smearing (mm)
     // TODO: update
     si_vtx_x = si_vtx_y = 600E-3;

     // pitch-induced error (mm), later adjusted by parameters.h
     // TODO: update
     si_de_P_L = si_de_P_R = 13E-3;

     // optics: v_x and v_y [1], L_x and L_y [mm]
     // sent by Frici on 23 Oct 2015
     v_x_L_2_F = -1.8975238180785; L_x_L_2_F = 0.10624114216534E3; v_y_L_2_F = -0.000000003186328; L_y_L_2_F = 261.86107319594E3;
     v_x_L_2_N = -2.1876926248587; L_x_L_2_N = 2.95354551535812E3; v_y_L_2_N = +0.020514691932280; L_y_L_2_N = 236.73917844622E3;
     v_x_L_1_F = -2.2756291135852; L_x_L_1_F = 3.81642926806849E3; v_y_L_1_F = +0.026731729097787; L_y_L_1_F = 229.12591622497E3;

     v_x_R_1_F = -2.2582096378676; L_x_R_1_F = 3.76173451557219E3; v_y_R_1_F = +0.026620752547344; L_y_R_1_F = 229.83172867404E3;
     v_x_R_2_N = -2.1682134167719; L_x_R_2_N = 2.89089335973313E3; v_y_R_2_N = +0.020429520698897; L_y_R_2_N = 237.53468452721E3;
     v_x_R_2_F = -1.8712479992497; L_x_R_2_F = 0.01733151135160E3; v_y_R_2_F = -0.000000023213780; L_y_R_2_F = 262.95254622452E3;


     // optics: x-y coupling
 /*
     la_x_L_F = la_x_L_N = la_x_R_N = la_x_R_F = 0.;    // mm
     la_y_L_F = la_y_L_N = la_y_R_N = la_y_R_F = 0.;    // mm
 */

     /*
     // optics imperfections
     double opt_cov_data[] = {
         1.66491785322919E-5,    7.89369350809322E-4,    -6.32104648991575E-5,    -2.59256651347955E-3,    1.32082198894547E-5,    6.74825862436010E-4,    -7.05099468507492E-5,    -2.90814857624182E-3,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,
         7.89369350809322E-4,    4.25168512900734E-2,    -2.69123774586626E-3,    -1.27194879518425E-1,    6.22063175217557E-4,    3.68812664207966E-2,    -3.00284426127789E-3,    -1.43008843891627E-1,   0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,
         -6.32104648991575E-5,    -2.69123774586626E-3,    3.87268155124952E-4,    1.16710801015928E-2,    -5.29597981974794E-5,    -2.04304942253293E-3,    4.31885407514716E-4,    1.29371751873752E-2,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,
         -2.59256651347955E-3,    -1.27194879518425E-1,    1.16710801015928E-2,    4.67495352905620E-1,    -2.06850163729353E-3,    -1.07850543630469E-1,    1.30191344764728E-2,    5.23596681087111E-1,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,
         1.32082198894547E-5,    6.22063175217557E-4,    -5.29597981974794E-5,    -2.06850163729353E-3,    1.05589617320310E-5,    5.24577953037806E-4,    -5.90732823172670E-5,    -2.31625003829467E-3,   0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,
         6.74825862436010E-4,    3.68812664207966E-2,    -2.04304942253293E-3,    -1.07850543630469E-1,    5.24577953037806E-4,    3.26401174262052E-2,    -2.27990513284437E-3,    -1.21628367533737E-1,   0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,
         -7.05099468507492E-5,    -3.00284426127789E-3,    4.31885407514716E-4,    1.30191344764728E-2,    -5.90732823172670E-5,    -2.27990513284437E-3,    4.81643176886755E-4,    1.44316029530475E-2,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,
         -2.90814857624182E-3,    -1.43008843891627E-1,    1.29371751873752E-2,    5.23596681087111E-1,    -2.31625003829467E-3,    -1.21628367533737E-1,    1.44316029530475E-2,    5.86636930463780E-1,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,

         0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,        1.70596E-005,    7.58302E-004,    -6.32105E-005,    -2.54045E-003,    1.37330E-005,    6.34320E-004,    -7.05009E-005,    -2.84011E-003,
         0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,        7.58302E-004,    4.53036E-002,    -2.69559E-003,    -1.31790E-001,    5.81606E-004,    4.04788E-002,    -3.00795E-003,    -1.48951E-001,
         0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,        -6.32105E-005,    -2.69559E-003,    3.87058E-004,    1.16688E-002,    -5.29702E-005,    -2.04924E-003,    4.31625E-004,    1.29347E-002,
         0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,        -2.54045E-003,    -1.31790E-001,    1.16688E-002,    4.74781E-001,    -2.00147E-003,    -1.13824E-001,    1.30165E-002,    5.33061E-001,
         0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,        1.37330E-005,    5.81606E-004,    -5.29702E-005,    -2.00147E-003,    1.12292E-005,    4.71804E-004,    -5.90750E-005,    -2.22875E-003,
         0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,        6.34320E-004,    4.04788E-002,    -2.04924E-003,    -1.13824E-001,    4.71804E-004,    3.72829E-002,    -2.28723E-003,    -1.29352E-001,
         0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,        -7.05009E-005,    -3.00795E-003,    4.31625E-004,    1.30165E-002,    -5.90750E-005,    -2.28723E-003,    4.81325E-004,    1.44289E-002,
         0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,    0.0000000000E+00,        -2.84011E-003,    -1.48951E-001,    1.29347E-002,    5.33061E-001,    -2.22875E-003,    -1.29352E-001,    1.44289E-002,    5.98935E-001
     };
     opt_cov.SetMatrixArray(opt_cov_data);

     TMatrixDSymEigen eig_decomp(opt_cov);
     TVectorD eig_values(eig_decomp.GetEigenValues());
     TMatrixDSym S(16);
     for (unsigned int i = 0; i < 16; i++)
         S(i, i) = (eig_values(i) >= 0.) ? sqrt(eig_values(i)) : 0.;
     opt_per_gen = eig_decomp.GetEigenVectors() * S;

     // alignment uncertainties
     si_de_x = 30E-3;
     si_de_y_R = 70E-3;
     si_de_y_D = 20E-3;
     si_tilt = 2E-3;

     // other uncertainties
     si_th_y_RL_assym_unc = 0.25;
     */
 }

 //----------------------------------------------------------------------------------------------------

 void Environment::UseMatchedOptics()
 {
 }

 /******************************************************************
  * common.h
  ******************************************************************/

  void Init(const std::string &dgnStr)
  {
      // ---------- apply standard settings ----------

      Init_base();

      diagonal = dUnknown;

      if (dgnStr.compare("45b_56t") == 0)
      {
          Init_45b_56t();
          diagonal = d45b_56t;
          th_y_sign = +1.;
      }

      if (dgnStr.compare("45t_56b") == 0)
      {
          Init_45t_56b();
          diagonal = d45t_56b;
          th_y_sign = -1.;
      }

      if (dgnStr.compare("combined") == 0)
      {
          diagonal = dCombined;
      }

      if (dgnStr.compare("45b_56b") == 0)
      {
          diagonal = ad45b_56b;

          Init_45t_56b();
          th_y_sign = -1.;
      }

      if (dgnStr.compare("45t_56t") == 0)
      {
          diagonal = ad45t_56t;

          Init_45b_56t();
          th_y_sign = +1.;
      }

      if (diagonal == dUnknown)
      {
          printf("Init > unrecognised diagonal %s\n", dgnStr.c_str());
      }

      // ---------- apply additional settings, if any ----------

  #ifdef USE_INIT_ADDITIONAL
      printf(">> Applying Init_additional\n");
      Init_additional();
  #endif

      // ---------- print important information ----------

      printf(">> bunchMap\n");
      // for (BunchMap::iterator rit = bunchMap.begin(); rit != bunchMap.end(); ++rit)
      // {
      //     printf("\trun %u: ", rit->first);
      //     for (unsigned int i = 0; i < rit->second.size(); i++)
      //     {
      //         if (i > 0)
      //             printf(", ");
      //         printf("%u", rit->second[i]);
      //     }
      //     printf("\n");
      // }
  }



 /******************************************************************
 * parameters.h
 ******************************************************************/

 void Init_base()
 {
     // load global settings
     Init_global();

     // directory for large-data storage (e.g. distilled ntuples)
     storageDir = "root://eostotem.cern.ch//eos/totem/user/j/jkaspar/analyses/elastic/6500GeV,beta90,10sigma/DS1/";

     // list of subdirectories with distilled ntuples
     distilledNtuples.push_back("block1");

     // selection of bunches
     keepAllBunches = true;
     //bunchMap[94882].push_back(0);

     // alignment settings
     /*
     AlignmentSource alSrc;
     alSrc.SetAlignmentA(atNone);
     alSrc.SetAlignmentB(atNone);
     alSrc.SetAlignmentC(atNone);

     alSrc.cnst.a_L_F = 0E-3; alSrc.cnst.b_L_F = 0E-3; alSrc.cnst.c_L_F = 0E-3;
     alSrc.cnst.a_L_N = 0E-3; alSrc.cnst.b_L_N = 0E-3; alSrc.cnst.c_L_N = 0E-3;
     alSrc.cnst.a_R_N = 0E-3; alSrc.cnst.b_R_N = 0E-3; alSrc.cnst.c_R_N = 0E-3;
     alSrc.cnst.a_R_F = 0E-3; alSrc.cnst.b_R_F = 0E-3; alSrc.cnst.c_R_F = 0E-3;

     alignmentSources.push_back(alSrc);
     */

 #if 0
     anal.use_time_dependent_resolutions = false;

     anal.use_3outof4_efficiency_fits = false;
     anal.use_pileup_efficiency_fits = false;
     anal.inefficiency_3outof4 = 0.;                // diagonal dependent
     anal.inefficiency_shower_near = 0.03;
     anal.inefficiency_pile_up = 0.;                // diagonal dependent
     anal.inefficiency_trigger = 0.;

     anal.bckg_corr = 1.;

     anal.L_int_eff = 0.;    // mb^-1, diagonal dependent

     anal.eff_th_y_min = 200E-6; // TODO

     anal.t_min_fit = 0.027; // TODO
 #endif

     anal.alignment_t0 = 0.;            // beginning of the first time-slice
     anal.alignment_ts = 5.*60.;        // time-slice in s

     anal.alignmentYRanges["L_2_F"] = Analysis::AlignmentYRange(-22.0, -10.0, 9.0, 22.0);
     anal.alignmentYRanges["L_2_N"] = Analysis::AlignmentYRange(-20.0, - 9.0, 8.4, 20.0);
     anal.alignmentYRanges["L_1_F"] = Analysis::AlignmentYRange(-19.3, - 8.8, 8.0, 19.6);

     anal.alignmentYRanges["R_1_F"] = Analysis::AlignmentYRange(-20.0, - 7.4, 8.4, 20.0);
     anal.alignmentYRanges["R_2_N"] = Analysis::AlignmentYRange(-20.5, - 7.6, 8.8, 20.8);
     anal.alignmentYRanges["R_2_F"] = Analysis::AlignmentYRange(-23.0, - 8.2, 9.6, 24.0);

 #if 0
     // TODO
     unsmearing_file = "";    // diagonal dependent
     //unsmearing_object = "cf,<binning>/exp3/corr_final";
     //unsmearing_object = "cf,<binning>/exp3+exp4/corr_final";
     unsmearing_object = "ff";

     // TODO
     luminosity_data_file = "../fill_3549_lumiCalc2.py_V04-02-04_lumibylsXing.csv";
 #endif
 }

 //----------------------------------------------------------------------------------------------------

 void Init_45b_56t()
 {
     Init_global_45b_56t();

     // analysis settings
     anal.cut1_a = 1.; anal.cut1_c = +0.5E-6; anal.cut1_si = 9.5E-6;
     anal.cut2_a = 1.; anal.cut2_c = -0.21E-6; anal.cut2_si = 2.8E-6;

     anal.cut5_a = 0.107200; anal.cut5_c = -0.010; anal.cut5_si = 0.016;
     anal.cut6_a = 0.105559; anal.cut6_c = -0.002; anal.cut6_si = 0.019;

     anal.cut7_a = 181.; anal.cut7_c = 0.; anal.cut7_si = 0.012;

 #if 0
     // TODO
     //unsmearing_file = "unfolding_fit_45b_56t_old.root";

     anal.inefficiency_3outof4 = 0.0; // TODO
     anal.inefficiency_pile_up = 0.0; // TODO

     anal.L_int_eff = 79.42E3;    // TODO
 #endif
 }

 //----------------------------------------------------------------------------------------------------

 void Init_45t_56b()
 {
     Init_global_45t_56b();

     // analysis settings
     anal.cut1_a = 1.; anal.cut1_c = +0.33E-6; anal.cut1_si = 10.0E-6;
     anal.cut2_a = 1.; anal.cut2_c = +0.35E-6; anal.cut2_si = 2.8E-6;

     anal.cut5_a = 0.10671; anal.cut5_c = 0.; anal.cut5_si = 0.018;
     anal.cut6_a = 0.10564; anal.cut6_c = 0.; anal.cut6_si = 0.018;

     anal.cut7_a = 179.; anal.cut7_c = 0.; anal.cut7_si = 0.012;

 #if 0
     // TODO
     //unsmearing_file = "unfolding_fit_45b_56t_old.root";

     anal.inefficiency_3outof4 = 0.0; // TODO
     anal.inefficiency_pile_up = 0.0; // TODO

     anal.L_int_eff = 79.42E3;    // TODO
 #endif
 }



 /******************************************************************
  * common_algorithms.h
  ******************************************************************/

 void BuildBinning(const Analysis &anal, const string &type, double* &binEdges, unsigned int &bins, bool verbose)
 {
     if (verbose)
         printf(">> BuildBinning(%s)\n", type.c_str());

     std::vector<double> be;
     double w;

     // same as in the low-|t| analysis
     if (type.compare("ub") == 0)
     {
         w = 2E-3;
         double t = 0.;
         while (t < anal.t_max_full)
         {
             be.push_back(t);
             t += w;
         }

         bins = be.size() - 1;
         binEdges = new double[bins + 1];
         for (unsigned int i = 0; i <= bins; i++)
             binEdges[i] = be[i];

         return;
     }

     // between t_min_full and t_min
     unsigned int N_bins_low = 4;
     w = (anal.t_min - anal.t_min_full) / N_bins_low;
     for (unsigned int i = 0; i < N_bins_low; i++)
         be.push_back(anal.t_min_full + w * i);

     // between t_min and t_max
     unsigned int N_bins_cen = 200;

     if (type.compare("eb") == 0)
     {
         double B = 3.;
         for (unsigned int bi = 0; bi < N_bins_cen; bi++)
             be.push_back( - log( (1. - double(bi) / N_bins_cen) * exp(-B*anal.t_min) + double(bi) * exp(-B*anal.t_max) / N_bins_cen ) / B );
     }

     if (type.find("ob") == 0)
     {
         // extract parameters
         size_t p1 = type.find("-", 0);
         size_t p2 = type.find("-", p1 + 1);
         size_t p3 = type.find("-", p2 + 1);

         double n_smearing_sigmas = atof(type.substr(p1+1, p2-p1-1).c_str());
         string stat_unc_label = type.substr(p2+1, p3-p2-1);
         double bs_max = atof(type.substr(p3+1).c_str());

         // load generators
         //TFile *f_in = TFile::Open("/afs/cern.ch/work/j/jkaspar/analyses/elastic/6500GeV,beta90,10sigma/binning/generators.root");
         TFile *f_in = TFile::Open("generators.root");
         TGraph *g_rms_t = (TGraph *) f_in->Get("g_rms_t");
         TGraph *g_bs_fsu = (TGraph *) f_in->Get( ("g_bs_stat_unc_" + stat_unc_label).c_str() );

         double t = anal.t_min;
         while (t < anal.t_max)
         {
             be.push_back(t);

             double w = max(n_smearing_sigmas * g_rms_t->Eval(t), g_bs_fsu->Eval(t));
             double t_c = t + w/2.;
             w = max(n_smearing_sigmas * g_rms_t->Eval(t_c), g_bs_fsu->Eval(t_c));
             if (w > bs_max)
                 w = bs_max;

             t += w;
         }

         delete f_in;
     }

     // between t_max and t_max_full
     unsigned int N_bins_high = 4;
     w = (anal.t_max_full - anal.t_max) / N_bins_high;
     for (unsigned int i = 0; i <= N_bins_high; i++)
         be.push_back(anal.t_max + w * i);

     // return results
     bins = be.size() - 1;
     binEdges = new double[be.size()];
     for (unsigned int i = 0; i < be.size(); i++)
     {
         binEdges[i] = be[i];
         if (verbose)
             printf("\tbi = %4u: %.4E\n", i, binEdges[i]);
     }
 }

 Binning BuildBinningRDF(const Analysis &anal, const string &type){
     Binning b;
     unsigned int N_bins;
     double *bin_edges;
     BuildBinning(anal, type, bin_edges, N_bins);
     b.N_bins = N_bins;
     b.bin_edges = bin_edges;
     return b;
 }

 //----------------------------------------------------------------------------------------------------

 bool CalculateAcceptanceCorrections(double th_y_sign,
         const Kinematics &k, const Analysis &anal,
         double &phi_corr, double &div_corr)
 {
     // ---------- smearing component ----------

     /*
     if ((k.th_x_L < anal.th_x_lcut_L) || (k.th_x_R < anal.th_x_lcut_R) || (k.th_x_L > anal.th_x_hcut_L) || (k.th_x_R > anal.th_x_hcut_R))
         return true;
     */

     if ((th_y_sign * k.th_y_L < anal.th_y_lcut_L) || (th_y_sign * k.th_y_R < anal.th_y_lcut_R)
         || (th_y_sign * k.th_y_L > anal.th_y_hcut_L) || (th_y_sign * k.th_y_R > anal.th_y_hcut_R)){
         //anal.Print();
         //printf("salida 1\n");
         return true;
     }

     /*
     double LB_x_L = anal.th_x_lcut_L - k.th_x, UB_x_L = anal.th_x_hcut_L - k.th_x;
     double LB_x_R = anal.th_x_lcut_R - k.th_x, UB_x_R = anal.th_x_hcut_R - k.th_x;
     double F_x_L = (UB_x_L > LB_x_L) ? ( TMath::Erf(UB_x_L / anal.si_th_x_1arm_L / sqrt(2.)) - TMath::Erf(LB_x_L / anal.si_th_x_1arm_L / sqrt(2.)) ) / 2. : 0.;
     double F_x_R = (UB_x_R > LB_x_R) ? ( TMath::Erf(UB_x_R / anal.si_th_x_1arm_R / sqrt(2.)) - TMath::Erf(LB_x_R / anal.si_th_x_1arm_R / sqrt(2.)) ) / 2. : 0.;
     double F_x = F_x_L * F_x_R;
     */
     double F_x = 1.;

     double th_y_abs = th_y_sign * k.th_y;

     double UB_y = min(anal.th_y_hcut_R - th_y_abs, th_y_abs - anal.th_y_lcut_L);
     double LB_y = max(anal.th_y_lcut_R - th_y_abs, th_y_abs - anal.th_y_hcut_L);
     double F_y = (UB_y > LB_y) ? ( TMath::Erf(UB_y / anal.si_th_y_1arm) - TMath::Erf(LB_y / anal.si_th_y_1arm) ) / 2. : 0.;

     //printf(">> F_x_L = %E, F_x_R = %E, F_y = %E\n", F_x_L, F_x_R, F_y);

     div_corr = 1./(F_x * F_y);

     // ---------- phi component ----------

     // apply safety margins to avoid excessive smearing component
     //double th_x_lcut = max(anal.th_x_lcut_L, anal.th_x_lcut_R) + 3.0E-6;
     //double th_x_hcut = min(anal.th_x_hcut_L, anal.th_x_hcut_R) - 3.0E-6;
     double th_x_lcut = anal.th_x_lcut;
     double th_x_hcut = anal.th_x_hcut;

     //double th_y_lcut = max(anal.th_y_lcut_L, anal.th_y_lcut_R) + 0.2E-6;
     //double th_y_hcut = min(anal.th_y_hcut_L, anal.th_y_hcut_R) - 1.0E-6;
     double th_y_lcut = anal.th_y_lcut;
     double th_y_hcut = anal.th_y_hcut;

     if (k.th_x <= th_x_lcut || k.th_x >= th_x_hcut || th_y_abs <= th_y_lcut || th_y_abs >= th_y_hcut){
         return true;
     }

     // get all intersections
     set<double> phis;

     if (k.th > th_y_lcut)
     {
         double phi = asin(th_y_lcut / k.th);
         double ta_x = k.th * cos(phi);
         if (th_x_lcut < ta_x && ta_x < th_x_hcut)
             phis.insert(phi);
         if (th_x_lcut < -ta_x && -ta_x < th_x_hcut)
             phis.insert(M_PI - phi);
     }

     if (k.th > th_y_hcut)
     {
         double phi = asin(th_y_hcut / k.th);
         double ta_x = k.th * cos(phi);
         if (th_x_lcut < ta_x && ta_x < th_x_hcut)
             phis.insert(phi);
         if (th_x_lcut < -ta_x && -ta_x < th_x_hcut)
             phis.insert(M_PI - phi);
     }

     if (k.th > fabs(th_x_hcut))
     {
         double phi = acos(fabs(th_x_hcut) / k.th);
         double ta_y = k.th * sin(phi);
         if (th_y_lcut < ta_y && ta_y < th_y_hcut)
             phis.insert(phi);
     }

     if (k.th > fabs(th_x_lcut))
     {
         double phi = acos(fabs(th_x_lcut) / k.th);
         double ta_y = k.th * sin(phi);
         if (th_y_lcut < ta_y && ta_y < th_y_hcut)
             phis.insert(M_PI - phi);
     }

     // the number of intersections must be even
     if ((phis.size() % 2) == 1)
     {
         printf("ERROR: odd number of intersections in acceptance calculation\n");
     }

     // no intersection => no acceptances
     if (phis.size() == 0){
         return true;
     }

     // calculate arc-length in within acceptance
     double phiSum = 0.;
     for (set<double>::iterator it = phis.begin(); it != phis.end(); ++it)
     {
         double phi_start = *it;
         ++it;
         double phi_end = *it;

         phiSum += phi_end - phi_start;
     }

     phi_corr = 2. * M_PI / phiSum;
     return false;
 }

 //----------------------------------------------------------------------------------------------------

 bool SkipRun(unsigned int /*run*/, unsigned int /*file*/, bool /*strict = true */)
 {
     return false;
 }

 //----------------------------------------------------------------------------------------------------

 // map: run number (8372) --> list of triggered bunches
 // typedef std::map<unsigned int, std::vector<unsigned int> > BunchMap;
 //
 // bool keepAllBunches;
 // BunchMap bunchMap;

 bool SkipBunch(unsigned int run, unsigned bunch)
 {
     if (keepAllBunches)
         return false;

     const std::vector<unsigned int> &bunches = bunchMap[run];

     return (find(bunches.begin(), bunches.end(), bunch) == bunches.end());
 }

 //----------------------------------------------------------------------------------------------------

 // returns the beam for which the bunch is non-colliding
 // for colliding bunches returns zero
 unsigned int NonCollidingBunch(unsigned int /*run*/, unsigned /*bunch*/)
 {
     /*
     if (run == 8318) {
         if (bunch == 994)
             return 1;
         if (bunch == 991)
             return 2;
     }

     if (run >= 8333 && run <= 8341)
     {
         if (bunch == 900)
             return 1;
         if (bunch == 991)
             return 2;
     }

     if (run >= 8367 && run <= 8372)
     {
         if (bunch == 3104 || bunch == 3130 || bunch == 3156 || bunch == 3078)
             return 1;
         if (bunch == 3143 || bunch == 3169 || bunch == 3195 || bunch == 3117)
             return 2;
     }
     */

     return 0;
 }

 //----------------------------------------------------------------------------------------------------

 bool IsZeroBias(unsigned int trigger, unsigned int /*run*/, unsigned int /*event*/)
 {
     return ((trigger & 512) != 0);
 }

 //----------------------------------------------------------------------------------------------------

 HitData ProtonTransport(const Kinematics & /*k*/, const Environment & /*env*/)
 {
     HitData h;

     // TODO
     /*
     h.x_L_F = -env.L_x_L_F*k.th_x_L + env.v_x_L_F*k.vtx_x   - env.la_x_L_F*k.th_y_L;
     h.y_L_F = -env.L_y_L_F*k.th_y_L + env.v_y_L_F*k.vtx_y   - env.la_y_L_F*k.th_x_L;

     h.x_L_N = -env.L_x_L_N*k.th_x_L + env.v_x_L_N*k.vtx_x   - env.la_x_L_N*k.th_y_L;
     h.y_L_N = -env.L_y_L_N*k.th_y_L + env.v_y_L_N*k.vtx_y   - env.la_y_L_N*k.th_x_L;

     h.x_R_N = +env.L_x_R_N*k.th_x_R + env.v_x_R_N*k.vtx_x   + env.la_x_R_N*k.th_y_R;
     h.y_R_N = +env.L_y_R_N*k.th_y_R + env.v_y_R_N*k.vtx_y   + env.la_y_R_N*k.th_x_R;

     h.x_R_F = +env.L_x_R_F*k.th_x_R + env.v_x_R_F*k.vtx_x   + env.la_x_R_F*k.th_y_R;
     h.y_R_F = +env.L_y_R_F*k.th_y_R + env.v_y_R_F*k.vtx_y   + env.la_y_R_F*k.th_x_R;
     */

     return h;
 }

 HitData ApplyFineAlignment( unsigned int &timestamp,
                             double &x_L_1_F, double &x_L_2_N, double &x_L_2_F,
                             double &x_R_1_F, double &x_R_2_N, double &x_R_2_F,
                             double &y_L_1_F, double &y_L_2_N, double &y_L_2_F,
                             double &y_R_1_F, double &y_R_2_N, double &y_R_2_F)
 {
     UnitHitData L_1_F, L_2_N, L_2_F;

     L_1_F.x = x_L_1_F; L_1_F.y = y_L_1_F; //L_1_F.x = x_L_1_F;
     L_2_N.x = x_L_2_N; L_2_N.y = y_L_2_N; // L_2_N
     L_2_F.x = x_L_2_F; L_2_F.y = y_L_2_F; // L_2_F

     UnitHitData R_1_F, R_2_N, R_2_F;

     R_1_F.x = x_R_1_F; R_1_F.y = y_R_1_F;
     R_2_N.x = x_R_2_N; R_2_N.y = y_R_2_N;
     R_2_F.x = x_R_2_F; R_2_F.y = y_R_2_F;

     HitData h_al;

     h_al.L_1_F = L_1_F;
     h_al.L_2_N = L_2_N;
     h_al.L_2_F = L_2_F;

     h_al.R_1_F = R_1_F;
     h_al.R_2_N = R_2_N;
     h_al.R_2_F = R_2_F;

     extern vector<AlignmentSource> alignmentSources;

     for (unsigned int i = 0; i < alignmentSources.size(); ++i)
     {
       AlignmentData alData = alignmentSources[i].Eval(timestamp);
       h_al = h_al.ApplyAlignment(alData);
     }

     return h_al;
 };

 Kinematics DoReconstruction(HitData &h)
 {
     Kinematics k;
     extern Environment env ;
     // single-arm kinematics reconstruction
     // th_x: linear regression
     // th_y: from hit positions
     // vtx_x: linear regression

     double D_x_L = - env.v_x_L_2_N * env.L_x_L_2_F + env.v_x_L_2_F * env.L_x_L_2_N;
     k.th_x_L = (env.v_x_L_2_N * h.L_2_F.x - env.v_x_L_2_F * h.L_2_N.x) / D_x_L;
     k.vtx_x_L = (- h.L_2_N.x * env.L_x_L_2_F + h.L_2_F.x * env.L_x_L_2_N) / D_x_L;

     double D_x_R = + env.v_x_R_2_N * env.L_x_R_2_F - env.v_x_R_2_F * env.L_x_R_2_N;
     k.th_x_R = (env.v_x_R_2_N * h.R_2_F.x - env.v_x_R_2_F * h.R_2_N.x) / D_x_R;
     k.vtx_x_R = (+ h.R_2_N.x * env.L_x_R_2_F - h.R_2_F.x * env.L_x_R_2_N) / D_x_R;

     double th_y_L_2_N = - h.L_2_N.y / env.L_y_L_2_N;
     double th_y_L_2_F = - h.L_2_F.y / env.L_y_L_2_F;
     k.th_y_L = (th_y_L_2_N + th_y_L_2_F) / 2.;

     double th_y_R_2_N = + h.R_2_N.y / env.L_y_R_2_N;
     double th_y_R_2_F = + h.R_2_F.y / env.L_y_R_2_F;
     k.th_y_R = (th_y_R_2_N + th_y_R_2_F) / 2.;

     double D_y_L = - env.v_y_L_2_N * env.L_y_L_2_F + env.v_y_L_2_F * env.L_y_L_2_N;
     //k.th_y_L = (env.v_y_L_2_N * L_2_F.y - env.v_y_L_2_F * L_2_N.y) / D_y_L;
     k.vtx_y_L = (- h.L_2_N.y * env.L_y_L_2_F + h.L_2_F.y * env.L_y_L_2_N) / D_y_L;

     double D_y_R = + env.v_y_R_2_N * env.L_y_R_2_F - env.v_y_R_2_F * env.L_y_R_2_N;
     //k.th_y_R = (env.v_y_R_2_N * R_2_F.y - env.v_y_R_2_F * R_2_N.y) / D_y_R;
     k.vtx_y_R = (+ h.R_2_N.y * env.L_y_R_2_F - h.R_2_F.y * env.L_y_R_2_N) / D_y_R;

     // double-arm kinematics reconstruction
     // th_x: from hit positions, L-R average
     // th_y: from hit positions, L-R average
     // vtx_x: from hit positions, L-R average

     k.th_x = (k.th_x_L + k.th_x_R) / 2.;
     k.th_y = (k.th_y_L + k.th_y_R) / 2.;

     k.vtx_x = (k.vtx_x_L + k.vtx_x_R) / 2.;
     k.vtx_y = (k.vtx_y_L + k.vtx_y_R) / 2.;

     // theta reconstruction
     double th_sq = k.th_x*k.th_x + k.th_y*k.th_y;
     k.th = sqrt(th_sq);
     k.phi = atan2(k.th_y, k.th_x);

     // t reconstruction
     k.t_x = env.p*env.p * k.th_x * k.th_x;
     k.t_y = env.p*env.p * k.th_y * k.th_y;
     k.t = k.t_x + k.t_y;

     return k;
 };


 Correction CalculateAcceptanceCorrectionsRDF(const Kinematics &k)
 {
     Correction correction;
     extern Analysis anal;
     extern double th_y_sign;
     double phi_corr = 0., div_corr = 0.;

     bool skip = CalculateAcceptanceCorrections(th_y_sign, k, anal, phi_corr, div_corr);

     correction.skip = skip;
     correction.phi_corr = phi_corr;
     correction.div_corr = div_corr;
     correction.corr = phi_corr * div_corr;
     return correction;
 };

 // Wrapper around anal.Skiptime
 bool SkipTime( unsigned int &timestamp){
     extern Analysis anal ;
     return anal.SkipTime(timestamp);
 };

 // Custom function to replace original check in line distributions.cc::820
 bool SkipTimeInterval( unsigned int &timestamp, int &tgd, int &tgr ){
     double time_group_interval = 1.;    // s
     int time_group = int(timestamp / time_group_interval);
     return  ( (time_group % tgd) != tgr);
 };

 // Custom function to replace original check in line distributions.cc::1021
 double getNorm_corr( unsigned int &timestamp ){
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
 };

 // Custom function to replace original check in line distributions.cc::1048
 double getNormalization( double &norm_corr ){
     extern Analysis anal;

     double normalization = anal.bckg_corr * norm_corr / anal.L_int;

     return normalization;
 };

 // FIXME Optimize this
 // This functions is meant to be used in a RDF::Define
 // where a column will be defined containing a 1 value for event
 double One(){
     return 1.;
 };

