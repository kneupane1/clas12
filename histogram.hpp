
TH2D *MomVsBeta =
    new TH2D("MomVsBeta", "Momentum vs Beta_all", 500, 0, 5, 500, 0, 1.2);

TH2D *MomVsBetaPostitve = new TH2D(
    "MomVsBetaPositive", "Momentum vs Beta_positive", 500, 0, 5, 500, 0, 1.2);

TH2D *MomVsBetaNegative = new TH2D(
    "MomVsBetaNegative", "Momentum vs Beta_negative", 500, 0, 5, 500, 0, 1.2);

TH2D *MomVsBetaProton = new TH2D("MomVsBetaProton", "Momentum vs Beta_Proton",
                                 500, 0, 5, 500, 0, 1.2);
TH2D *MomVsBetaPip =
    new TH2D("MomVsBetaPip", "Momentum vs Beta_Pip", 500, 0, 5, 500, 0, 1.2);
TH2D *MomVsBetaPim =
    new TH2D("MomVsBetaPim", "Momentum vs Beta_Pim", 500, 0, 5, 500, 0, 1.2);

TH2D *wq2 = new TH2D("wq2", "W vs Q^{2}", 500, 0, 5, 500, 0, 6.0);
TH1D *w = new TH1D("w", "W", 500, 0, 5);
TH1D *mm_pim = new TH1D("mm_pim", "mm_pim", 500, -4.9, 4.9);

TH2D *sf_hist =
    new TH2D("sf_hist", "Electron Sampling Fraction", 500, 0, 11, 500, 0, 0.5);

TH1D *PCAL_hist = new TH1D(" PCAL", "PCAL ; PCAL (GeV); Entries", 500, 0, 3.6);

TH1D *ECAL_hist = new TH1D(" ECAL", "ECAL ; ECAL (GeV); Entries", 500, 0, 2.6);
TH2D *pcal_SF_vs_P_hist =
    new TH2D("pcal_SF_vs_P", "pcal_SF_vs_P; P (GeV);	pcal_SF ", 500, 0.0, 11,
             500, 0.0, 0.55);
TH2D *ecal_SF_vs_P_hist =
    new TH2D("ecal_SF_vs_P", "ecal_SF_vs_P; P (GeV);	ecal_SF ", 500, 0.0, 11,
             500, 0.0, 0.55);
TH2D *ECAL_vs_PCAL_hist =
    new TH2D("ECAL_vs_PCAL", "ECAL_vs_PCAL; PCAL (GeV);	ECAL (GeV)", 500, 0,
             1.5, 500, 0, 1.6);
TH2D *ECout_vs_ECin_hist =
    new TH2D("ECout_vs_ECin", "ECout_vs_ECin; ECin (GeV);	ECout (Gev) ",
             500, 0, 1.2, 500, 0, 1);
TH2D *pcal_Fiducial_cut_hist =
    new TH2D("pcal_Fiducial_cut", "pcal_Fiducial_cut; x/cm;	y/cm ", 1500,
             -400, 400, 1500, -400, 400);

TH2D *deltaT_1a_prot =
    new TH2D("deltaT_1a_prot", "#Deltat Proton", 500, 0, 7.0, 500, -10, 10);
TH2D *deltaT_1a_pip =
    new TH2D("deltaT_1a_pion", "#Deltat #pi^{+}", 500, 0, 7.0, 500, -10, 10);
TH2D *deltaT_1a_pim =
    new TH2D("deltaT_1a_pion_m", "#Deltat #pi^{-}", 500, 0, 7.0, 500, -10, 10);

TH2D *deltaT_1b_prot =
    new TH2D("deltaT_1b_prot", "#Deltat Proton", 500, 0, 7.0, 500, -10, 10);
TH2D *deltaT_1b_pip =
    new TH2D("deltaT_1b_pion", "#Deltat #pi^{+}", 500, 0, 7.0, 500, -10, 10);
TH2D *deltaT_1b_pim =
    new TH2D("deltaT_1b_pion_m", "#Deltat #pi^{-}", 500, 0, 7.0, 500, -10, 10);

TH2D *deltaT_2_prot =
    new TH2D("deltaT_2_prot", "#Deltat Proton", 500, 0, 7.0, 500, -10, 10);
TH2D *deltaT_2_pip =
    new TH2D("deltaT_2_pion", "#Deltat #pi^{+}", 500, 0, 7.0, 500, -10, 10);
TH2D *deltaT_2_pim =
    new TH2D("deltaT_2_pion_m", "#Deltat #pi^{-}", 500, 0, 7.0, 500, -10, 10);

TH2D *deltaT_ctof_prot =
    new TH2D("deltaT_ctof_prot", "#Deltat Proton", 500, 0, 7.0, 500, -10, 10);
TH2D *deltaT_ctof_pip =
    new TH2D("deltaT_ctof_pion", "#Deltat #pi^{+}", 500, 0, 7.0, 500, -10, 10);
TH2D *deltaT_ctof_pim = new TH2D("deltaT_ctof_pion_m", "#Deltat #pi^{-}", 500,
                                 0, 7.0, 500, -10, 10);
