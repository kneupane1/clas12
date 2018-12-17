#include "delta_t.hpp"
#
double Q2_calc(TLorentzVector e_mu, TLorentzVector e_mu_prime) {
  TLorentzVector q_mu = (e_mu - e_mu_prime);
  return -q_mu.Mag2();
}
//	Calcualting W
//	Gotten from s channel [(gamma - P)^2 == s == w^2]
//	Sqrt√[M_p^2 - Q^2 + 2 M_p gamma]
double W_calc(TLorentzVector e_mu, TLorentzVector e_mu_prime) {
  TLorentzVector q_mu = (e_mu - e_mu_prime);
  TVector3 p_mu_3(0, 0, 0);
  TLorentzVector p_mu;
  p_mu.SetVectM(p_mu_3, MASS_P);
  return (p_mu + q_mu).Mag();
}
double missing_mass_calc(TLorentzVector e_mu_prime, TLorentzVector p_mu_prime,
                         TLorentzVector pip_mu_prime,
                         TLorentzVector pim_mu_prime) {

  TLorentzVector e_mu(0.0, 0.0, 10.6, 10.6);
  TVector3 p_mu_3(0, 0, 0);
  TLorentzVector p_mu;
  //  TLorentzVector p_mu_prime;
  p_mu.SetVectM(p_mu_3, MASS_P);
  // std::cout << "p_mu" << p_mu.Mag() << '\n';

  return (e_mu + p_mu - e_mu_prime - p_mu_prime -
          pip_mu_prime /*- pip_mu_prime*/)
      .Mag();
}
