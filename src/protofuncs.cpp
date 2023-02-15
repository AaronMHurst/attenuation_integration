#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "constants.hpp"
#include "protolib.hpp"

double Prototypes::ReturnValue()
{
    return *RESULT;
}

double Prototypes::getMuNeutron(double A, double& rho, double elemental_sigmaNABS, double temperature, double E, double v, double reduced_lambda, bool zsolt_method)
{
    double result = 0;

    if (zsolt_method == true) {
        //Check Zsolt's spreadsheet; I think there is a mistake in the expression
        double this_part_is_ok = (1 / A) * 0.6023 * elemental_sigmaNABS * rho;
        double this_part_is_different = sqrt(T_th / temperature);
        result = this_part_is_ok * this_part_is_different;
    }

    else {
        //I think this is the correct routine below

        double B = T_th / temperature;
        //double B = sqrt(E_th / E);
        //double B = v_th / v;
        //double B = reduced_lambda;

        //All four above methods should give identical evaluation of B factor

        result = (Avagadro_NA / A) * rho * Convert_b_cm2 * elemental_sigmaNABS * B;
    }

    RESULT = &result;
    return result;
}

double Prototypes::getMuGamma(double& rho, double u_rho)
{
    double result = u_rho * rho;
    RESULT = &result;
    return result;
}

double Prototypes::calcAttenuation(double thickness, double u_gamma, double u_neutron, double alpha_rad, int cmORmm, bool uGammaOnly)
{

    double K = 0.0;
    if (uGammaOnly == false) {
        K = (u_gamma / cos(alpha_rad)) + (u_neutron / sin(alpha_rad));
        //std::cout << 1/K <<std::endl;
    } else if (uGammaOnly == true) {
        //std::cout<<"hello "<<uGammaOnly<<std::endl;
        K = (u_gamma / cos(alpha_rad));
        //std::cout << 1/K <<std::endl;
    }

    double function_exp = 0;
    double ratio = 0;
    double result = 0;

    switch (cmORmm) {
    case 1:
        //Use method below for [mm]

        function_exp = exp(-(thickness / 10.0) * K);
        ratio = (1 - function_exp) / K;
        result = ratio / (thickness / 10.0);
        break;

    case 2:
        //Use method below for [cm]

        function_exp = exp(-thickness * K);
        ratio = (1 - function_exp) / K;
        result = ratio / thickness;
        break;

    default:
        break;
    }

    RESULT = &result;
    return result;
}

double Prototypes::getErrMuGamma(double& rho, double u_rho)
{
    double u_gamma = u_rho * rho;
    double d_u_gamma = (u_gamma / 100) * 5.0;

    double result = d_u_gamma;
    RESULT = &result;
    return result;
}

double Prototypes::getErrAttenuation(double thickness, double alpha_rad, double I_I0, int cmORmm, double u_neutron, double u_gamma, double d_u_neutron, double d_u_gamma, bool uGammaOnly)
{
    double result = 0;

    switch (cmORmm) {
    case 1:
        thickness = thickness / 10.0;
        break;
    case 2:
        thickness = thickness;
        break;
    default:
        break;
    }

    double sum_u_coeffs = 0.0;
    if (uGammaOnly == false)
        sum_u_coeffs = u_gamma / cos(alpha_rad) + u_neutron / sin(alpha_rad);
    else if (uGammaOnly == true)
        sum_u_coeffs = u_gamma / cos(alpha_rad);

    double func_exp = thickness * exp((-thickness) * sum_u_coeffs);

    double common_factor = (1.0 / sum_u_coeffs) * (func_exp - I_I0);
    common_factor = pow(common_factor, 2.0);

    double sqsum_d_u = 0.0;
    if (uGammaOnly == false) {
        d_u_neutron = d_u_neutron / sin(alpha_rad);
        d_u_gamma = d_u_gamma / cos(alpha_rad);
        sqsum_d_u = pow(d_u_neutron, 2.0) + pow(d_u_gamma, 2.0);
    }

    else if (uGammaOnly == true) {
        d_u_gamma = d_u_gamma / cos(alpha_rad);
        sqsum_d_u = pow(d_u_gamma, 2.0);
    }

    result = sqrt(common_factor * sqsum_d_u);

    RESULT = &result;
    return result;
}

//The following functions are used in the analysis of compound samples
double Prototypes::getRMM(int num_atoms, double A)
{
    double result = num_atoms * A;
    RESULT = &result;
    return result;
}

double Prototypes::getWeightedMF(double u_rho, double w_mass_factor)
{
    double result = u_rho * w_mass_factor;
    RESULT = &result;
    return result;
}

double Prototypes::getMuNeutronCompound(double rho_compound, double temperature, double ratio_CS_A, double w_mass_factor, bool zsolt_method)
{
    double constants = Avagadro_NA * Convert_b_cm2;
    double input_vars = rho_compound;
    if (zsolt_method == true)
        input_vars = input_vars * sqrt(T_th / temperature);
    else
        input_vars = input_vars * (T_th / temperature);
    double loop_vars = ratio_CS_A * w_mass_factor;
    double result = constants * input_vars * loop_vars;
    RESULT = &result;
    return result;
}

//Runtime prompts

void Prototypes::printZsoltMethod()
{
    std::cout << std::endl;
    //std::cout << std::setfill('#') << std::setw(45) << std::endl;
    std::cout << "################################################" << std::endl;
    std::cout << "# The following option is for debugging only!  #" << std::endl;
    std::cout << "# Generally you should answer \'No\' (i.e. 0)    #" << std::endl;
    std::cout << "################################################" << std::endl;
    std::cout << std::endl;

    return;
}

void Prototypes::printAbsorber(std::string absorber)
{
    std::cout << std::endl;
    std::cout << "################################################" << std::endl;
    std::cout << "# Calculating attenuation coefficients for: "
              << absorber << std::endl;
    std::cout << "################################################" << std::endl;
    std::cout << std::endl;

    return;
}

void Prototypes::printNISTCalculator()
{
    std::cout << std::endl;
    std::cout << "##########################################################" << std::endl;
    std::cout << "# Above result neglects contribution due to incoherent   #" << std::endl;
    std::cout << "# neutron scattering. Compare above neutron coefficient  #" << std::endl;
    std::cout << "# to NIST calculator:                                    #" << std::endl;
    std::cout << "# http://www.ncnr.nist.gov/instruments/bt1/neutron.html  #" << std::endl;
    std::cout << "##########################################################" << std::endl;

    return;
}
