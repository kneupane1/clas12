#include "TCanvas.h"
#include "TChain.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "delta_t.hpp"
#include "functions.hpp"
#include "histogram.hpp"
#include <TFile.h>
#include <chrono>
#include <cstdlib>
#include <iostream>

std::vector<int> *pid;
std::vector<float> *p;
std::vector<float> *px;
std::vector<float> *py;
std::vector<float> *pz;
std::vector<float> *beta;
std::vector<int> *charge;
std::vector<int> *ec_pcal_sec;

std::vector<float> *ec_tot_energy;
std::vector<float> *ec_pcal_energy;
std::vector<float> *ec_ecin_energy;
std::vector<float> *ec_ecout_energy;
std::vector<float> *ec_pcal_x;
std::vector<float> *ec_pcal_y;

std::vector<float> *sc_ftof_1a_time;
std::vector<float> *sc_ftof_1a_path;

std::vector<float> *sc_ftof_1b_time;
std::vector<float> *sc_ftof_1b_path;

std::vector<float> *sc_ftof_2_time;
std::vector<float> *sc_ftof_2_path;

std::vector<float> *sc_ctof_time;
std::vector<float> *sc_ctof_path;
static const double Pival = TMath::Pi();

// Calcuating Q^2
// q^mu^2 = (e^mu - e^mu')^2 = -Q^2
TLorentzVector e_mu(0.0, 0.0, BEAM, BEAM);
TLorentzVector e_mu_prime(0.0, 0.0, 0.0, 0.0);
TLorentzVector p_mu_prime(0.0, 0.0, 0.0, 0.0);
TLorentzVector pip_mu_prime(0.0, 0.0, 0.0, 0.0);
TLorentzVector pim_mu_prime(0.0, 0.0, 0.0, 0.0);
double mm = 0, w_cal = 0;
double energy = 0, PCAL = 0, ECAL = 0, ECin = 0, ECout = 0, pcal_SF = 0,
       ecal_SF = 0;
TChain *clas12 = new TChain("clas12", "clas12");

void test(char *fin, char *fout) {
  TFile *out = new TFile(fout, "RECREATE");

  clas12->Add(fin);
  clas12->SetBranchAddress("pid", &pid);
  clas12->SetBranchAddress("p", &p);
  clas12->SetBranchAddress("px", &px);
  clas12->SetBranchAddress("py", &py);
  clas12->SetBranchAddress("pz", &pz);
  clas12->SetBranchAddress("beta", &beta);
  clas12->SetBranchAddress("charge", &charge);
  clas12->SetBranchAddress("ec_tot_energy", &ec_tot_energy);
  clas12->SetBranchAddress("ec_pcal_sec", &ec_pcal_sec);
  clas12->SetBranchAddress("ec_pcal_x", &ec_pcal_x);
  clas12->SetBranchAddress("ec_pcal_y", &ec_pcal_y);
  clas12->SetBranchAddress("ec_pcal_energy", &ec_pcal_energy);
  clas12->SetBranchAddress("ec_ecin_energy", &ec_ecin_energy);
  clas12->SetBranchAddress("ec_ecout_energy", &ec_ecout_energy);

  clas12->SetBranchAddress("sc_ftof_1a_time", &sc_ftof_1a_time);
  clas12->SetBranchAddress("sc_ftof_1a_path", &sc_ftof_1a_path);
  clas12->SetBranchAddress("sc_ftof_1b_time", &sc_ftof_1b_time);
  clas12->SetBranchAddress("sc_ftof_1b_path", &sc_ftof_1b_path);
  clas12->SetBranchAddress("sc_ftof_2_time", &sc_ftof_2_time);
  clas12->SetBranchAddress("sc_ftof_2_path", &sc_ftof_2_path);

  clas12->SetBranchAddress("sc_ctof_time", &sc_ctof_time);
  clas12->SetBranchAddress("sc_ctof_path", &sc_ctof_path);

  int num_of_events = (int)clas12->GetEntries();
  auto start_full = std::chrono::high_resolution_clock::now();
  for (int current_event = 0; current_event < num_of_events; current_event++) {
    clas12->GetEntry(current_event);
    if (pid->size() == 0 || pid->at(0) != 11)
      continue;
    if (pid->at(0) != 22 && pid->at(0) != 0 && p->at(0) != 0) {

      int sec_PCAL = ec_pcal_sec->at(0) - 1;
      double x_PCAL = ec_pcal_x->at(0);
      double y_PCAL = ec_pcal_y->at(0);
      double x_PCAL_rot = y_PCAL * sin(sec_PCAL * 60.0 * Pival / 180) +
                          x_PCAL * cos(sec_PCAL * 60.0 * Pival / 180);
      double y_PCAL_rot = y_PCAL * cos(sec_PCAL * 60.0 * Pival / 180) -
                          x_PCAL * sin(sec_PCAL * 60.0 * Pival / 180);
      double angle_PCAL = 60;
      double height_PCAL = 45;
      double slope_PCAL = 1 / tan(0.5 * angle_PCAL * Pival / 180);
      double left_PCAL = (height_PCAL - slope_PCAL * y_PCAL_rot);
      double right_PCAL = (height_PCAL + slope_PCAL * y_PCAL_rot);
      double radius2_PCAL = pow(height_PCAL + 6, 2) - pow(y_PCAL_rot, 2);
      if (x_PCAL_rot > left_PCAL && x_PCAL_rot > right_PCAL &&
          pow(x_PCAL_rot, 2) > radius2_PCAL && x_PCAL_rot < 362) {
        pcal_Fiducial_cut_hist->Fill(x_PCAL, y_PCAL);
      }
      e_mu_prime.SetXYZM(px->at(0), py->at(0), pz->at(0), MASS_E);
      //  std::cout << "e_mu_prime 1111  " << e_mu_prime.P() << '\n';

      energy = ec_tot_energy->at(0);
      PCAL = ec_pcal_energy->at(0);
      ECin = ec_ecin_energy->at(0);
      ECout = ec_ecout_energy->at(0);
      ECAL = ECin + ECout;
      if (ec_tot_energy->at(0) != 0) {
        sf_hist->Fill(p->at(0), ec_tot_energy->at(0) / p->at(0));
      }
      if (PCAL != 0) {
        PCAL_hist->Fill(PCAL);
      }
      if (ECAL != 0) {
        ECAL_hist->Fill(ECAL);
      }
      if (ECAL != 0) {
        if (PCAL > 0.01)
          ECAL_vs_PCAL_hist->Fill(PCAL, ECAL);
      }
      if (ECout != 0) {
        if (ECin > 0.01)
          ECout_vs_ECin_hist->Fill(ECin, ECout);
      }
      pcal_SF = PCAL / e_mu_prime.P();
      if (pcal_SF != 0) {
        pcal_SF_vs_P_hist->Fill(e_mu_prime.P(), pcal_SF);
      }
      ecal_SF = ECAL / e_mu_prime.P();
      //  if (ecal_SF > 0.05 && ecal_SF < 0.8) {
      if (ecal_SF != 0) {
        ecal_SF_vs_P_hist->Fill(e_mu_prime.P(), ecal_SF);
      }
    }

    for (size_t i = 1; i < pid->size(); i++) {
      if (beta->at(i) < 0.05 || charge->at(i) == 0)
        continue;
      MomVsBeta->Fill(p->at(i), beta->at(i));

      if (charge->at(i) > 0) {
        MomVsBetaPostitve->Fill(p->at(i), beta->at(i));
        if (pid->at(i) == 2212) {
          MomVsBetaProton->Fill(p->at(i), beta->at(i));
        } else if (pid->at(i) == 211) {
          MomVsBetaPip->Fill(p->at(i), beta->at(i));
        }

      } else if (charge->at(i) < 0) {
        MomVsBetaNegative->Fill(p->at(i), beta->at(i));
        if (pid->at(i) == -211) {
          MomVsBetaPim->Fill(p->at(i), beta->at(i));
        }
      }
    }

    double sf = (ec_tot_energy->at(0) / p->at(0));
    if (sf < 0.3 && sf > 0.2) {
      e_mu_prime.SetXYZM(px->at(0), py->at(0), pz->at(0), MASS_E);
      w_cal = W_calc(e_mu, e_mu_prime);
      wq2->Fill(W_calc(e_mu, e_mu_prime), Q2_calc(e_mu, e_mu_prime));
      w->Fill(W_calc(e_mu, e_mu_prime));
    }

    double vertex = 0.0;
    double beta = 0.0;

    //  std::cout << "p_mu_prime_11111 " << p_mu_prime.P() << '\n';

    if (sc_ftof_1a_time->at(0) == sc_ftof_1a_time->at(0)) {
      vertex = vertex_time(sc_ftof_1a_time->at(0), sc_ftof_1a_path->at(0), 1.0);
    } else if (sc_ftof_1b_time->at(0) == sc_ftof_1b_time->at(0)) {
      vertex = vertex_time(sc_ftof_1b_time->at(0), sc_ftof_1b_path->at(0), 1.0);
    } else {
      continue;
    }
    for (size_t part = 0; part < pid->size(); part++) {
      if (p->at(part) == 0)
        continue;

      //  continue;
      if (charge->at(part) == 1) {
        beta =
            1.0 / sqrt(1.0 + (MASS_P / p->at(part)) * (MASS_P / p->at(part)));
        dt_P = vertex - vertex_time(sc_ftof_1a_time->at(part),
                                    sc_ftof_1a_path->at(part), beta);
        if (dt_P == dt_P) {

          deltaT_1a_prot->Fill(p->at(part), dt_P);
        }
        beta = 1.0 /
               sqrt(1.0 + (MASS_PIP / p->at(part)) * (MASS_PIP / p->at(part)));
        dt_PIP = vertex - vertex_time(sc_ftof_1a_time->at(part),
                                      sc_ftof_1a_path->at(part), beta);
        if (dt_PIP == dt_PIP) {

          deltaT_1a_pip->Fill(p->at(part), dt_PIP);
        }

        beta =
            1.0 / sqrt(1.0 + (MASS_P / p->at(part)) * (MASS_P / p->at(part)));
        dt_P = vertex - vertex_time(sc_ftof_1b_time->at(part),
                                    sc_ftof_1b_path->at(part), beta);
        if (dt_P == dt_P) {
          deltaT_1b_prot->Fill(p->at(part), dt_P);
        }
        beta = 1.0 /
               sqrt(1.0 + (MASS_PIP / p->at(part)) * (MASS_PIP / p->at(part)));
        dt_PIP = vertex - vertex_time(sc_ftof_1b_time->at(part),
                                      sc_ftof_1b_path->at(part), beta);
        if (dt_PIP == dt_PIP) {
          deltaT_1b_pip->Fill(p->at(part), dt_PIP);
        }

        beta =
            1.0 / sqrt(1.0 + (MASS_P / p->at(part)) * (MASS_P / p->at(part)));
        dt_P = vertex - vertex_time(sc_ftof_2_time->at(part),
                                    sc_ftof_2_path->at(part), beta);
        if (dt_P == dt_P) {
          deltaT_2_prot->Fill(p->at(part), dt_P);
        }

        beta = 1.0 /
               sqrt(1.0 + (MASS_PIP / p->at(part)) * (MASS_PIP / p->at(part)));
        dt_PIP = vertex - vertex_time(sc_ftof_2_time->at(part),
                                      sc_ftof_2_path->at(part), beta);
        if (dt_PIP == dt_PIP) {
          deltaT_2_pip->Fill(p->at(part), dt_PIP);
        }

        //  {
        beta =
            1.0 / sqrt(1.0 + (MASS_P / p->at(part)) * (MASS_P / p->at(part)));
        dt_P = vertex - vertex_time(sc_ctof_time->at(part),
                                    sc_ctof_path->at(part), beta);
        if (dt_P == dt_P && dt_P > -1.0 && dt_P < 0.5) {
          deltaT_ctof_prot->Fill(p->at(part), dt_P);
          p_mu_prime.SetXYZM(px->at(part), py->at(part), pz->at(part), MASS_P);
        }

        beta = 1.0 /
               sqrt(1.0 + (MASS_PIP / p->at(part)) * (MASS_PIP / p->at(part)));
        dt_PIP = vertex - vertex_time(sc_ctof_time->at(part),
                                      sc_ctof_path->at(part), beta);
        if (dt_PIP == dt_PIP && dt_PIP < 0.6 && dt_PIP > -0.5) {
          deltaT_ctof_pip->Fill(p->at(part), dt_PIP);
          pip_mu_prime.SetXYZM(px->at(part), py->at(part), pz->at(part),
                               MASS_PIP);
        }

        //}
      } else if (charge->at(part) == -1) {

        beta = 1.0 /
               sqrt(1.0 + (MASS_PIP / p->at(part)) * (MASS_PIP / p->at(part)));
        dt_PIP = vertex - vertex_time(sc_ftof_1a_time->at(part),
                                      sc_ftof_1a_path->at(part), beta);
        if (dt_PIP == dt_PIP)
          deltaT_1a_pim->Fill(p->at(part), dt_PIP);

        beta = 1.0 /
               sqrt(1.0 + (MASS_PIP / p->at(part)) * (MASS_PIP / p->at(part)));
        dt_PIP = vertex - vertex_time(sc_ftof_1b_time->at(part),
                                      sc_ftof_1b_path->at(part), beta);
        if (dt_PIP == dt_PIP)
          deltaT_1b_pim->Fill(p->at(part), dt_PIP);

        beta = 1.0 /
               sqrt(1.0 + (MASS_PIP / p->at(part)) * (MASS_PIP / p->at(part)));
        dt_PIP = vertex - vertex_time(sc_ftof_2_time->at(part),
                                      sc_ftof_2_path->at(part), beta);
        if (dt_PIP == dt_PIP)
          deltaT_2_pim->Fill(p->at(part), dt_PIP);

        beta = 1.0 /
               sqrt(1.0 + (MASS_PIP / p->at(part)) * (MASS_PIP / p->at(part)));
        dt_PIP = vertex - vertex_time(sc_ctof_time->at(part),
                                      sc_ctof_path->at(part), beta);
        double sf = (ec_tot_energy->at(0) / p->at(0));
        if (sf < 0.3 && sf > 0.2) {
          e_mu_prime.SetXYZM(px->at(0), py->at(0), pz->at(0), MASS_E);
          w_cal = W_calc(e_mu, e_mu_prime);
          if (w_cal > 1.2 && w_cal < 1.4) {
            if (dt_PIP == dt_PIP) {
              deltaT_ctof_pim->Fill(p->at(part), dt_PIP);
              pim_mu_prime.SetXYZM(px->at(part), py->at(part), pz->at(part),
                                   MASS_PIP);
            }
          }
        }
      }
    }
    if (p_mu_prime.P() != 0 && pip_mu_prime.P() != 0 && w_cal > 1.2 &&
        w_cal < 1.4) {
      mm =
          missing_mass_calc(e_mu_prime, p_mu_prime, pip_mu_prime, pim_mu_prime);

      mm_pim->Fill(mm);
    }
  }

  out->cd();
  TDirectory *wvsq2 = out->mkdir("folder");
  TCanvas *c1 = new TCanvas("WvsQ2", "WvsQ2");
  wvsq2->cd();
  c1->cd();
  wq2->Draw("colz");
  c1->Write();

  TCanvas *c2 = new TCanvas("w", "w", 10, 10, 900, 600);
  c2->Divide(2, 1);
  c2->cd(1);
  w->Draw("same");
  c2->cd(2);
  mm_pim->Draw("same");
  c2->Write();

  TCanvas *c3 = new TCanvas("dt_ftof", "dt_ftof", 10, 10, 900, 600);
  c3->Divide(3, 3);
  c3->cd(1);
  deltaT_1a_prot->Draw("colz");
  c3->cd(2);
  deltaT_1a_pip->Draw("colz");
  c3->cd(3);
  deltaT_1a_pim->Draw("colz");
  c3->cd(4);
  deltaT_1b_prot->Draw("colz");
  c3->cd(5);
  deltaT_1b_pip->Draw("colz");
  c3->cd(6);
  deltaT_1b_pim->Draw("colz");
  c3->cd(7);
  deltaT_2_prot->Draw("colz");
  c3->cd(8);
  deltaT_2_pip->Draw("colz");
  c3->cd(9);
  deltaT_2_pim->Draw("colz");
  c3->Write();

  TCanvas *c4 = new TCanvas("dt_ctof", "dt_ctof", 910, 10, 900, 600);
  c4->Divide(3);
  c4->cd(1);
  deltaT_ctof_prot->Draw("colz");
  c4->cd(2);
  deltaT_ctof_pip->Draw("colz");
  c4->cd(3);
  deltaT_ctof_pim->Draw("colz");
  c4->Write();

  TCanvas *c5 = new TCanvas("MomVsBeta", "MomVsBeta");
  c5->Divide(3, 2);
  c5->cd(1);
  MomVsBeta->Draw("colz");
  c5->cd(2);
  MomVsBetaPostitve->Draw("colz");
  c5->cd(3);
  MomVsBetaNegative->Draw("colz");
  c5->cd(4);
  MomVsBetaProton->Draw("colz");
  c5->cd(5);
  MomVsBetaPip->Draw("colz");
  c5->cd(6);
  MomVsBetaPim->Draw("colz");
  c5->Write();

  TDirectory *sf = out->mkdir("sf");
  sf->cd();
  sf_hist->SetOption("COLZ");
  sf_hist->Write();
  PCAL_hist->Write();
  ECAL_hist->Write();
  pcal_SF_vs_P_hist->SetOption("COLZ");
  pcal_SF_vs_P_hist->Write();
  ecal_SF_vs_P_hist->SetOption("COLZ");
  ecal_SF_vs_P_hist->Write();
  ECAL_vs_PCAL_hist->SetOption("COLZ");
  ECAL_vs_PCAL_hist->Write();
  ECout_vs_ECin_hist->SetOption("COLZ");
  ECout_vs_ECin_hist->Write();
  pcal_Fiducial_cut_hist->SetOption("COLZ");
  pcal_Fiducial_cut_hist->Write();
}
