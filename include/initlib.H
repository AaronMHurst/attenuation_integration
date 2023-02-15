#ifndef INITLIB_H
#define INITLIB_H

class Initialize{
public:
  Initialize(int num_elements, int num_atoms, double rho_compound, double temperature, double thickness, double user_cs, double user_err_cs, double user_A);
  Initialize();
  ~Initialize();

  void setNumberElements(int num_elements);
  void setNumberAtoms(int num_atoms, char* absorber);
  void setCompoundDensity(double rho_compound);
  void setTemperature(double temperature);
  void setThickness(double thickness);
  void setUserCSAndErrAndA(double user_cs, double user_err_cs, double user_A);

  int getNumberElements() const;
  int getNumberAtoms() const;
  double getCompoundDensity() const;
  double getTemperature() const;
  double getThickness() const;
  double getUserCS();
  double getUserErrCS();
  double getUserA();

  float getLambdaLambda0();
  float getWavelength();
  float getBeamVelocity();
  float getBeamEnergy();
  float calcDeg2Rad() const;

private:
  int NUMBER_OF_ELEMENTS;
  int NUMBER_OF_ATOMS;
  double RhoCompound;
  double Temperature;
  double Thickness;
  float ReducedWavelength;
  double UserCrossSection;
  double UserErrCrossSection;
  double UserA;
};


class Flag{
public:
  Flag(int cmORmm, int choice, int USE_ADOPTED_SIGMA, bool element);
  Flag();
  ~Flag();

  void chooseCMorMM(int cmORmm);
  void chooseElementOrCompound(int choice);
  void chooseAdoptedSigma(int USE_ADOPTED_SIGMA);

  int applyCMorMM() const;
  int applyElementOrCompound() const;
  int applyAdoptedSigma() const;

  bool elementTrueOrFalse(bool element);
  bool uGammaOnlyTrueOrFalse();

private:
  int CM_or_MM;
  int Choice;
  int UseAdoptedSigma;
  bool Element;
};

#endif
