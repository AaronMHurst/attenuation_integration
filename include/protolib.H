#ifndef PROTOLIB_H
#define PROTOLIB_H

class Prototypes{
public:
  Prototypes() {}
  ~Prototypes() {}

  double getMuNeutron(double A, double &rho, double elemental_sigmaNABS, double temperature, double E, double v, double reduced_lambda, bool zsolt_method);
  double getMuGamma(double &rho, double u_rho);
  double calcAttenuation(double thickness, double u_gamma, double u_neutron, double alpha_rad, int cmORmm, bool uGammaOnly);
  double getErrMuGamma(double &rho, double u_rho);
  double getErrAttenuation(double thickness, double alpha_rad, double I_I0, int cmORmm, double u_neutron, double u_gamma, double d_u_neutron, double d_u_gamma, bool uGammaOnly);
  double getRMM(int num_atoms, double A);
  double getWeightedMF(double u_rho, double w_mass_factor);
  double getMuNeutronCompound(double rho_compound, double temperature, double ratio_CS_A, double w_mass_factor, bool zsolt_method);

  double ReturnValue();

  void printZsoltMethod();
  void printNISTCalculator();
  void printAbsorber(std::string absorber);

private:
  double *RESULT;
};

#endif
