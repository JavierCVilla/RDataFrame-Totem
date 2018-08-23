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

