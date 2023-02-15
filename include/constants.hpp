#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

// Budapest Research Reactor parameters
const double alpha_deg = 30.0;

// Thermal-neutron kinematics
const double T_th = 293.0;
const double v_th = 2200.0;
const double E_th = 0.02526;
const double lambda_0 = 1.8;

// Fundamental constants
// Not all are being used!

const double PI = 3.141592654;
const double k_SI = 1.3807E-23; // J/K
const double k_nuclear = 8.617E-05; // eV/K

const double h_SI = 6.626E-34; // J.s
const double h_nuclear = 4.136E-15; // eV.s

const double Avagadro_NA = 6.0221367E+23; // mol^{-1}

const double Convert_b_cm2 = 1.0E-24; // 1b = 1 x 10^{-24} cm^2

// Some array definitions
#define MAX_NUMBER_ELEMENTS 10
#define MAX_NUMBER_LINES 500

#endif
