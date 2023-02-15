#include <cmath>
#include <cstdlib>
#include <iostream>

#include "constants.hpp"
#include "initlib.hpp"

Initialize::Initialize(int num_elements, int num_atoms, double rho_compound, double temperature, double thickness, double user_cs, double user_err_cs, double user_A)
{
    NUMBER_OF_ELEMENTS = num_elements;
    NUMBER_OF_ATOMS = num_atoms;
    RhoCompound = rho_compound;
    Temperature = temperature;
    Thickness = thickness;
    UserCrossSection = user_cs;
    UserErrCrossSection = user_err_cs;
    UserA = user_A;
}

Initialize::Initialize()
{
    NUMBER_OF_ELEMENTS = 0;
    NUMBER_OF_ATOMS = 0;
    RhoCompound = 0;
    Temperature = 0;
    Thickness = 0;
    UserCrossSection = 0;
    UserErrCrossSection = 0;
    UserA = 0;
}

Initialize::~Initialize()
{
}

void Initialize::setNumberElements(int num_elements)
{
    std::cout << "Number of Elements in compound?" << std::endl;
    std::cin >> num_elements;

    NUMBER_OF_ELEMENTS = num_elements;

    return;
}

int Initialize::getNumberElements() const
{
    return NUMBER_OF_ELEMENTS;
}

void Initialize::setNumberAtoms(int num_atoms, char* absorber)
{
    std::cout << "Number of atoms belonging to " << absorber << " (i.e. stoichiometry) ?" << std::endl;
    std::cin >> num_atoms;

    NUMBER_OF_ATOMS = num_atoms;

    return;
}

int Initialize::getNumberAtoms() const
{
    return NUMBER_OF_ATOMS;
}

void Initialize::setCompoundDensity(double rho_compound)
{
    std::cout << "Density of compound [g/cm^{3}]?" << std::endl;
    std::cin >> rho_compound;

    RhoCompound = rho_compound;

    return;
}

double Initialize::getCompoundDensity() const
{
    return RhoCompound;
}

void Initialize::setTemperature(double temperature)
{
    std::cout << "Temperature of neutron beam [K] ?" << std::endl;
    std::cin >> temperature;

    Temperature = temperature;

    return;
}

double Initialize::getTemperature() const
{
    return Temperature;
}

void Initialize::setThickness(double thickness)
{
    std::cin >> thickness;

    Thickness = thickness;

    return;
}

double Initialize::getThickness() const
{
    return Thickness;
}

void Initialize::setUserCSAndErrAndA(double user_cs, double user_err_cs, double user_A)
{
    std::cout << "Enter absorption cross section [b]: " << std::endl;
    std::cin >> user_cs;
    std::cout << "Enter uncertainty [b]: " << std::endl;
    std::cin >> user_err_cs;
    std::cout << "Enter mass (A) of isotope: " << std::endl;
    std::cin >> user_A;

    UserCrossSection = user_cs;
    UserErrCrossSection = user_err_cs;
    UserA = user_A;

    return;
}

double Initialize::getUserCS()
{
    return UserCrossSection;
}

double Initialize::getUserErrCS()
{
    return UserErrCrossSection;
}

double Initialize::getUserA()
{
    return UserA;
}

float Initialize::getLambdaLambda0()
{
    double result = T_th / Temperature;
    ReducedWavelength = result;
    return result;
}

float Initialize::getWavelength()
{
    double result = ReducedWavelength * lambda_0;
    return result;
}

float Initialize::getBeamVelocity()
{
    double result = v_th / ReducedWavelength;
    return result;
}

float Initialize::getBeamEnergy()
{
    double result = E_th / pow(ReducedWavelength, 2.0);
    return result;
}

float Initialize::calcDeg2Rad() const
{
    double result = (alpha_deg / 180.0) * PI;
    return result;
}

Flag::Flag(int cmORmm, int choice, int USE_ADOPTED_SIGMA, bool element)
{
    CM_or_MM = cmORmm;
    Choice = choice;
    UseAdoptedSigma = USE_ADOPTED_SIGMA;
    Element = element;
}

Flag::Flag()
{
    CM_or_MM = 0;
    Choice = 0;
    UseAdoptedSigma = 0;
    Element = 0;
}

Flag::~Flag()
{
}

void Flag::chooseCMorMM(int cmORmm)
{
    std::cout << "Sample in [mm] or [cm] ?\n1 - [mm]\n2 - [cm]" << std::endl;
    std::cin >> cmORmm;

    CM_or_MM = cmORmm;

    return;
}

int Flag::applyCMorMM() const
{
    return CM_or_MM;
}

void Flag::chooseElementOrCompound(int choice)
{
    std::cout << "Natural Element (1) or Compound Sample (2) ?" << std::endl;
    std::cin >> choice;

    Choice = choice;

    return;
}

int Flag::applyElementOrCompound() const
{
    return Choice;
}

bool Flag::elementTrueOrFalse(bool element)
{
    if (Choice == 1)
        element = true;
    else if (Choice == 2)
        element = false;
    else {
        std::cerr << "Invalid choice \nAbort" << std::endl;
        exit(-1);
    }

    return element;
}

void Flag::chooseAdoptedSigma(int USE_ADOPTED_SIGMA)
{
    std::cout << "Use adopted elemental absorption cross section from Mughabghab's\nAtlas of Neutron Resonances (Ed. 2006)?" << std::endl
              << std::endl;
    std::cout << "[Answer \'No\' if using an enriched isotope and provide isotopic-\nabsorption cross section and mass when prompted.]" << std::endl;
    std::cout << "1 - Yes\n2 - No" << std::endl;
    std::cin >> USE_ADOPTED_SIGMA;

    UseAdoptedSigma = USE_ADOPTED_SIGMA;

    return;
}

int Flag::applyAdoptedSigma() const
{
    return UseAdoptedSigma;
}

bool Flag::uGammaOnlyTrueOrFalse()
{
    int attenuation_type = 0;
    bool uGammaOnly = false;

    std::cout << "Calculate attenuation assuming coefficients for: " << std::endl;
    std::cout << "1 - gamma-ray attenuation only" << std::endl;
    std::cout << "2 - gamma-ray and neutron attenuation combined" << std::endl;

    std::cin >> attenuation_type;

    if (attenuation_type == 1)
        uGammaOnly = true;
    else if (attenuation_type == 2)
        uGammaOnly = false;
    else {
        std::cerr << "Invalid choice \nAbort" << std::endl;
        exit(-1);
    }

    return uGammaOnly;
}
