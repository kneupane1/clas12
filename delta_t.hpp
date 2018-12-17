#ifndef DT_H_GUARD
#define DT_H_GUARD

double BEAM = 10.6;
static const double MASS_P = 0.93827203;
static const double MASS_N = 0.93956556;
static const double MASS_E = 0.000511;
static const double MASS_PIP = 0.13957018;
static const double MASS_PIM = 0.13957018;
static const double MASS_PI0 = 0.1349766;
static const double MASS_KP = 0.493677;
static const double MASS_KM = 0.493677;
static const double MASS_G = 0.0;
static const double MASS_OMEGA = 0.78265;

const double c_special_units = 29.9792458;

double dt_P, dt_PIP, dt_E;

double vertex_time(double sc_time, double sc_pathlength,
                   double relatavistic_beta) {
  return sc_time - sc_pathlength / (relatavistic_beta * c_special_units);
}
#endif
