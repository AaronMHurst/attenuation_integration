//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//Program to calculate self-attenuation factors according to the linear     //
//combination of both gamma-ray- and neutron-attenuation coefficients       //
//integrated over the sample thickness.  Code generates attenuation factors //
//(I/I0) as a function of gamma-ray energy for pure (natural elemental)     //
//samples or stoichiometric compounds.                                      //
//                                                                          //
//The neutron-attenuation coefficient is determined using an assumed value  //
//for the neutron-capture cross section; the contribution from incoherent-  //
//neutron scattering (usually much smaller) is neglected in this            //
//calculation. These results can be compared to values obtained from the    //
//NIST pages available at:                                                  //
//http://www.ncnr.nist.gov/instruments/bt1/neutron.html                     //
//                                                                          //
//The overall uncertainty on the calculated attenuation factors (I/I0)      //
//is determined by assuming the conservative estimate of 5-% uncertainty    //
//on the gamma-ray mass-attenuation factors that are used to generate the   //
//energy-dependent gamma-ray attenuation coefficients, and appropriately    //
//propagating through this uncertainty together with the uncertainty on     //
//the neutron attenuation coefficient.                                      //
//                                                                          //
//Author: A. M. Hurst, LBNL                                                 //
//Date: June 2013                                                           //
//Email: AMHurst@lbl.gov                                                    //
//Update: May 15, 2015                                                      //
//Separated [u_g] and [u_g + u_n] options for attenuation calculation       //
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "constants.H"
#include "initlib.H"
#include "inputlib.H"
#include "protolib.H"

using namespace std;

int main()
{
    //Initialize member data of the Flag class using its constructor
    Flag(0, 0, 0, 0);

    //Instanciating a 'Decision'-type object of the Flag class
    Flag Decision;

    Initialize(0, 0, 0, 0, 0, 0, 0, 0);
    Initialize Number;

    Prototypes();
    Prototypes Function;

    //Declare u_gamma or [u_gamma + u_neutron] calculation
    bool uGammaOnly = Decision.uGammaOnlyTrueOrFalse();

    int cmORmm = 0;
    Decision.chooseCMorMM(cmORmm);
    cmORmm = Decision.applyCMorMM();

    double thickness = 0;

    switch (cmORmm) {
    case 1:
        cout << "Give sample thickness [mm]: " << endl;
        Number.setThickness(thickness);
        thickness = Number.getThickness();
        break;
    case 2:
        cout << "Give sample thickness [cm]: " << endl;
        Number.setThickness(thickness);
        thickness = Number.getThickness();
        break;
    default:
        cerr << "Invalid choice!" << endl;
        cerr << "\nABORT CALCULATION" << endl;
        exit(1);
        break;
    }

    double temperature = 0;
    Number.setTemperature(temperature);
    temperature = Number.getTemperature();

    float reduced_lambda = Number.getLambdaLambda0();
    cout << "Reduced neutron wavelength = " << reduced_lambda << endl;

    float lambda = Number.getWavelength();
    cout << "Neutron wavelength = " << lambda << " A" << endl;

    float beam_velocity = Number.getBeamVelocity();
    cout << "Neutron beam velocity = " << beam_velocity << " m/s" << endl;

    float beam_energy = Number.getBeamEnergy();
    cout << "Neutron beam energy = " << beam_energy << " eV" << endl;

    float alpha_rad = Number.calcDeg2Rad();
    cout << "Angle (alpha) of sample relative to beam direction = " << alpha_deg << " deg. = " << alpha_rad << " rad." << endl;

    ifstream infile2;
    infile2.open("phys_data/physical_data.dat");

    string phys_data_file = "physical_data.dat";

    std::vector<physical_data> in2;
    physical_data buffer2;

    int mlines = 0;

    while (true) {

        infile2 >> buffer2.chemical_symbol >> buffer2.atomic_number >> buffer2.A >> buffer2.rho >> buffer2.elemental_sigmaNABS >> buffer2.err_elemental_sigmaNABS >> buffer2.incoh_scatter_cs >> buffer2.err_incoh_scatter_cs;

        in2.push_back(buffer2);

        if (infile2.eof() || infile2.fail())
            break;

        mlines++;
    }
    in2.pop_back();
    cout << endl;
    cout << "Number of lines read from file " << phys_data_file << ": " << mlines << endl;
    cout << endl;
    infile2.close();

    ofstream outfile;
    outfile.open("self_attenuation.dat");

    outfile << "Neutron Beam Kinematics:" << endl
            << endl;
    outfile << "Angle (alpha) of sample relative to beam direction = " << alpha_deg << " deg. = " << alpha_rad << " rad." << endl;
    outfile << "Neutron beam temperature: " << temperature << " K" << endl;
    outfile << "Reduced neutron wavelength: " << reduced_lambda << endl;
    outfile << "Neutron wavelength: " << lambda << " A" << endl;
    outfile << "Neutron beam velocity: " << beam_velocity << " m/s" << endl;
    outfile << "Neutron beam energy: " << beam_energy << " eV" << endl;
    outfile << endl;

    int choice = 0;
    Decision.chooseElementOrCompound(choice);
    choice = Decision.applyElementOrCompound();

    //////////////////////////////////////////////////////////////////////
    //Temporary code snippet: compare calculation of the linear neutron //
    //attenuation coefficient with Zsolt's spreadsheet                  //
    //////////////////////////////////////////////////////////////////////

    bool zsolt_method = false;
    //Function.printZsoltMethod();
    //cout << "Use formalism in Zsolt's spreadsheet for calculation of \nlinear neutron-attenuation coefficient?\n0 - No \n1 - Yes" << endl;
    //cin >> zsolt_method;

    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////

    bool element = false;
    element = Decision.elementTrueOrFalse(element);

    if (element == true) {
        char filename[256];
        string path_to_file = "mass_data/";
        string end_filename = "_mu_rho.dat";
        string begin_filename;
        cout << "Chemical symbol for absorber?" << endl;
        cin >> begin_filename;

        Function.printAbsorber(begin_filename);

        string file = path_to_file + begin_filename + end_filename;

        sprintf(filename, file.c_str(), 256);

        ifstream infile;
        infile.open(filename);

        if (!infile) {
            cerr << filename << " does not exist! \nMake sure the case lettering is correct: \ne.g. write \'W\' rather than \'w\', or \n\'Pb\' rather than \'pb\' " << endl;
            cerr << "\nABORT CALCULATION" << endl;
            exit(1);
        } else
            cout << "Gamma-ray absorption coefficients for " << begin_filename << " taken from file:\n"
                 << filename << endl;

        string element_symbol = begin_filename;

        std::vector<record> in;
        record buffer;

        int nlines = 0;

        while (true) {

            infile >> buffer.energy >> buffer.u_rho;

            in.push_back(buffer);

            if (infile.eof() || infile.fail())
                break;

            nlines++;
        }
        in.pop_back();
        cout << "Number of lines read: " << nlines << endl
             << endl;
        infile.close();

        int* Z = 0;
        double A = 0, *rho = 0;
        double elemental_sigmaNABS = 0, err_elemental_sigmaNABS = 0;

        int USE_ADOPTED_SIGMA = 0;
        Decision.chooseAdoptedSigma(USE_ADOPTED_SIGMA);
        USE_ADOPTED_SIGMA = Decision.applyAdoptedSigma();

        bool FIND_CHEMICAL_SYMBOL = false;
        switch (USE_ADOPTED_SIGMA) {
        case 1:
            for (std::vector<physical_data>::iterator j = in2.begin(); j != in2.end(); ++j) {
                if (element_symbol == j->chemical_symbol) {
                    Z = &j->atomic_number;
                    A = j->A;
                    rho = &j->rho;
                    elemental_sigmaNABS = j->elemental_sigmaNABS;
                    err_elemental_sigmaNABS = j->err_elemental_sigmaNABS;

                    FIND_CHEMICAL_SYMBOL = true;
                }
            }
            break;

        case 2:
            for (std::vector<physical_data>::iterator j = in2.begin(); j != in2.end(); ++j) {
                if (element_symbol == j->chemical_symbol) {
                    Z = &j->atomic_number;
                    rho = &j->rho;

                    FIND_CHEMICAL_SYMBOL = true;
                }
            }

            Number.setUserCSAndErrAndA(elemental_sigmaNABS, err_elemental_sigmaNABS, A);
            elemental_sigmaNABS = Number.getUserCS();
            err_elemental_sigmaNABS = Number.getUserErrCS();
            A = Number.getUserA();

            break;

        default:
            FIND_CHEMICAL_SYMBOL = false;
            break;
        }

        if (FIND_CHEMICAL_SYMBOL == false) {
            cerr << "\nWrong choice!" << endl;
            cerr << "Element: \"" << element_symbol << "\" does not exist!" << endl;
            cerr << "Make sure the case lettering is correct: \ne.g. write \'W\' rather than \'w\', or \n\'Pb\' rather than \'pb\' " << endl;
            cerr << "\nABORT CALCULATION" << endl;
            exit(1);
        }

        outfile << "Absorber properties: " << endl;
        outfile << "Irradiated absorber: " << element_symbol << endl;
        outfile << "Element (Z) = " << *Z << endl;
        outfile << "Relative Atomic Number (A) = " << A << endl;

        switch (cmORmm) {
        case 1:
            outfile << "Sample thickness = " << thickness << " mm" << endl;
            break;
        case 2:
            outfile << "Sample thickness = " << thickness << " cm" << endl;
            break;
        default:
            break;
        }

        outfile << "Density = " << *rho << " g/cm^3" << endl;
        outfile << "Elemental neutron absorption cross section = "
                << elemental_sigmaNABS << " +/- " << err_elemental_sigmaNABS
                << " b" << endl;

        double u_neutron = 0.0, d_u_neutron = 0.0;
        if (uGammaOnly == true) {
            u_neutron = 0.0;
            d_u_neutron = 0.0;
        }

        else if (uGammaOnly == false) {
            u_neutron = Function.getMuNeutron(A, *rho, elemental_sigmaNABS, temperature, beam_energy, beam_velocity, reduced_lambda, zsolt_method);

            d_u_neutron = Function.getMuNeutron(A, *rho, err_elemental_sigmaNABS, temperature, beam_energy, beam_velocity, reduced_lambda, zsolt_method);
        }

        if (uGammaOnly == true) {
            outfile << endl
                    << "Neutron attenuation coefficient: N/A " << endl;

            cout << endl
                 << "Neutron attenuation coefficient: N/A " << endl;
        }

        else if (uGammaOnly == false) {
            outfile << endl
                    << "Neutron attenuation coefficient: " << u_neutron << " +/- " << d_u_neutron << " cm^{-1}" << endl;

            cout << endl
                 << "Neutron attenuation coefficient: " << u_neutron << " +/- " << d_u_neutron << " cm^{-1}" << endl;
            Function.printNISTCalculator();
        }

        //////////////////////////////////////////////////////////////////////
        //Temporary code snippet: compare calculation of the linear neutron //
        //attenuation coefficient with Zsolt's spreadsheet                  //
        //////////////////////////////////////////////////////////////////////
        if (zsolt_method == true) {
            outfile << "ZSOLT'S ROUTINE USED IN ABOVE CALCULATION" << endl;
            double u_neutron_rho_zsolt = u_neutron / (*rho);

            outfile << "cf. spreadsheet quantity: u_n/rho = "
                    << u_neutron_rho_zsolt << " cm^2/g" << endl;
        }
        //////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////

        outfile << endl;

        outfile << setiosflags(ios::left);
        outfile << setw(15) << "Energy"
                << "\t" << setw(20) << "u_g [cm^{-1}]"
                << "\t" << setw(20) << "I/I0"
                << "\t" << setw(20) << "dI_I0" << endl;
        outfile << endl;

        int LINE_COUNT = 0;
        for (std::vector<record>::iterator i = in.begin(); i != in.end(); i++) {
            LINE_COUNT++;
            int energy(i->energy);
            double u_rho(i->u_rho);

            double u_gamma = Function.getMuGamma(*rho, u_rho);

            double I_I0 = Function.calcAttenuation(thickness, u_gamma, u_neutron, alpha_rad, cmORmm, uGammaOnly);

            double d_u_gamma = Function.getErrMuGamma(*rho, u_rho);

            double d_I_I0 = Function.getErrAttenuation(thickness, alpha_rad, I_I0, cmORmm, u_neutron, u_gamma, d_u_neutron, d_u_gamma, uGammaOnly);

            //Nice output-formatting settings
            outfile << setiosflags(ios::left);
            //outfile<<setiosflags(ios::fixed);
            outfile << setiosflags(ios::scientific);
            outfile << setw(15) << energy << "\t" << setw(20) << setprecision(9) << u_gamma << "\t" << I_I0 << "\t" << setw(20) << d_I_I0 << endl;
        }

    }

    else if (element == false) {
        double rho_compound = 0;
        double compound_mu_gamma[MAX_NUMBER_LINES] = { 0 };
        //int energy_list[MAX_NUMBER_LINES] = {0};
        int* energy_list[MAX_NUMBER_LINES] = { 0 };
        double compound_mu_neutron = 0;

        Number.setCompoundDensity(rho_compound);
        rho_compound = Number.getCompoundDensity();

        int num_elements = 0;
        int num_atoms = 0;

        Number.setNumberElements(num_elements);
        num_elements = Number.getNumberElements();

        double RMM = 0;
        double RAM[MAX_NUMBER_ELEMENTS] = { 0 };

        int Z_LIST[MAX_NUMBER_ELEMENTS] = { 0 };
        double RHO_LIST[MAX_NUMBER_ELEMENTS] = { 0 };
        double CS_LIST[MAX_NUMBER_ELEMENTS] = { 0 };
        double d_CS_LIST[MAX_NUMBER_ELEMENTS] = { 0 };

        double ratio_CS_A[MAX_NUMBER_ELEMENTS] = { 0 };
        double err_ratio_CS_A[MAX_NUMBER_ELEMENTS] = { 0 };

        stringstream ss_absorber;
        stringstream ss_element[MAX_NUMBER_ELEMENTS];
        string str_element[MAX_NUMBER_ELEMENTS];
        string str_pure_element[MAX_NUMBER_ELEMENTS];

        int ELEMENT_COUNTER = 0;
        for (int i = 0; i < num_elements; i++) {
            ELEMENT_COUNTER++;

            string begin_filename;
            cout << "Chemical symbol for Absorber No. " << ELEMENT_COUNTER
                 << " ?" << endl;
            cin >> begin_filename;

            str_pure_element[i] = begin_filename;
            str_element[i] = begin_filename;
            ss_element[i] << str_element[i];

            char absorber[8];
            string element_symbol = begin_filename;
            sprintf(absorber, element_symbol.c_str(), 8);

            ////////////////////////////////////////////////////
            //Determining Mass-Weighting factors for absorbers//
            ////////////////////////////////////////////////////

            Number.setNumberAtoms(num_atoms, absorber);
            num_atoms = Number.getNumberAtoms();

            ss_absorber << absorber;
            if (num_atoms > 1) {
                ss_absorber << num_atoms;
                ss_element[i] << num_atoms;
                str_element[i] = ss_element[i].str();
            }

            int Z = 0;
            double A = 0, rho = 0;
            double elemental_sigmaNABS = 0, err_elemental_sigmaNABS = 0;

            int USE_ADOPTED_SIGMA = 0;
            Decision.chooseAdoptedSigma(USE_ADOPTED_SIGMA);
            USE_ADOPTED_SIGMA = Decision.applyAdoptedSigma();

            bool FIND_CHEMICAL_SYMBOL = false;
            switch (USE_ADOPTED_SIGMA) {
            case 1:
                for (std::vector<physical_data>::iterator j = in2.begin(); j != in2.end(); ++j) {
                    if (element_symbol == j->chemical_symbol) {
                        Z = j->atomic_number;
                        A = j->A;
                        rho = j->rho;
                        elemental_sigmaNABS = j->elemental_sigmaNABS;
                        err_elemental_sigmaNABS = j->err_elemental_sigmaNABS;

                        Z_LIST[i] = Z;
                        RHO_LIST[i] = rho;
                        CS_LIST[i] = elemental_sigmaNABS;
                        d_CS_LIST[i] = err_elemental_sigmaNABS;

                        ratio_CS_A[i] = (elemental_sigmaNABS / A);
                        err_ratio_CS_A[i] = (err_elemental_sigmaNABS / A);

                        FIND_CHEMICAL_SYMBOL = true;
                    }
                }
                break;

            case 2:
                for (std::vector<physical_data>::iterator j = in2.begin(); j != in2.end(); ++j) {
                    if (element_symbol == j->chemical_symbol) {
                        Z = j->atomic_number;
                        rho = j->rho;

                        Z_LIST[i] = Z;
                        RHO_LIST[i] = rho;

                        FIND_CHEMICAL_SYMBOL = true;
                    }
                }

                Number.setUserCSAndErrAndA(elemental_sigmaNABS, err_elemental_sigmaNABS, A);
                elemental_sigmaNABS = Number.getUserCS();
                err_elemental_sigmaNABS = Number.getUserErrCS();
                A = Number.getUserA();

                CS_LIST[i] = elemental_sigmaNABS;
                d_CS_LIST[i] = err_elemental_sigmaNABS;

                ratio_CS_A[i] = (elemental_sigmaNABS / A);
                err_ratio_CS_A[i] = (err_elemental_sigmaNABS / A);

                break;

            default:
                FIND_CHEMICAL_SYMBOL = false;
                break;
            }

            if (FIND_CHEMICAL_SYMBOL == false) {
                cerr << "\nWrong choice!" << endl;
                cerr << "Element: \"" << element_symbol << "\" does not exist!" << endl;
                cerr << "Make sure the case lettering is correct: \ne.g. write \'W\' rather than \'w\', or \n\'Pb\' rather than \'pb\' " << endl;
                cerr << "\nABORT CALCULATION" << endl;
                exit(1);
            }

            RMM += Function.getRMM(num_atoms, A);
            RAM[i] = Function.getRMM(num_atoms, A);
        }

        string str_absorber = ss_absorber.str();

        thickness = Number.getThickness();

        outfile << "Absorber properties: " << endl;
        outfile << "Irradiated compound absorber: " << str_absorber << endl;

        Function.printAbsorber(str_absorber);

        switch (cmORmm) {
        case 1:
            outfile << "Sample thickness = " << thickness << " mm" << endl;
            break;
        case 2:
            outfile << "Sample thickness = " << thickness << " cm" << endl;
            break;
        default:
            break;
        }

        cout << "RMM (" << str_absorber << ") = " << RMM << endl;
        outfile << "RMM (" << str_absorber << ") = " << RMM << endl;
        for (int i = 0; i < num_elements; i++) {
            cout << "RAM (" << str_element[i] << ") = " << RAM[i] << endl;
            outfile << "RAM (" << str_element[i] << ") = " << RAM[i] << endl;
            outfile << "Z (" << str_pure_element[i] << ") = " << Z_LIST[i] << endl;
        }

        double w_mass_factor = 0;
        double err_compound_mu_neutron = 0;

        for (int i = 0; i < num_elements; i++) {
            w_mass_factor = RAM[i] / RMM;

            outfile << "Weighted mass factor for " << str_element[i]
                    << " = " << w_mass_factor << endl;

            compound_mu_neutron += Function.getMuNeutronCompound(rho_compound, temperature, ratio_CS_A[i], w_mass_factor, zsolt_method);

            err_compound_mu_neutron += pow(Function.getMuNeutronCompound(rho_compound, temperature, err_ratio_CS_A[i], w_mass_factor, zsolt_method), 2.0);
        }
        err_compound_mu_neutron = sqrt(err_compound_mu_neutron);

        //Physical data lists here
        for (int i = 0; i < num_elements; i++) {
            outfile << "Density of " << str_pure_element[i] << " = " << RHO_LIST[i] << " g/cm^{3}" << endl;
            outfile << "Elemental neutron absorption cross section for " << str_pure_element[i] << " = " << CS_LIST[i] << " +/- " << d_CS_LIST[i] << " b" << endl;
        }

        outfile << endl
                << "Density of " << str_absorber << " compound used to generate linear attenuation \ncoefficients: " << Number.getCompoundDensity() << " g/cm^{3}" << endl;

        cout << endl
             << "Compound neutron attenuation coefficient: " << compound_mu_neutron << " +/- " << err_compound_mu_neutron << " cm^{-1}" << endl;
        Function.printNISTCalculator();
        cout << endl;

        outfile << endl
                << "Compound neutron attenuation coefficient: " << compound_mu_neutron << " +/- " << err_compound_mu_neutron << " cm^{-1}" << endl;
        //////////////////////////////////////////////////////////////////////
        //Temporary code snippet: compare calculation of the linear neutron //
        //attenuation coefficient with Zsolt's spreadsheet                  //
        //////////////////////////////////////////////////////////////////////
        if (zsolt_method == true) {
            outfile << "ZSOLT'S ROUTINE USED IN ABOVE CALCULATION" << endl;
            double u_neutron_rho_zsolt = compound_mu_neutron / (rho_compound);

            outfile << "cf. spreadsheet quantity: u_n/rho = "
                    << u_neutron_rho_zsolt << " cm^2/g" << endl;
        }
        //////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////
        outfile << endl;

        int nlines = 0;
        ELEMENT_COUNTER = 0;
        for (int i = 0; i < num_elements; i++) {
            ELEMENT_COUNTER++;

            w_mass_factor = RAM[i] / RMM;

            char filename[256];
            string path_to_file = "mass_data/";
            string end_filename = "_mu_rho.dat";
            string begin_filename;

            for (std::vector<physical_data>::iterator j = in2.begin(); j != in2.end(); ++j)
                if (Z_LIST[i] == j->atomic_number)
                    begin_filename = j->chemical_symbol;

            string file = path_to_file + begin_filename + end_filename;
            sprintf(filename, file.c_str(), 256);

            ifstream infile;
            infile.open(filename);

            std::vector<record> in;
            record buffer;

            if (!infile) {
                cerr << filename << " does not exist! \nMake sure the case lettering is correct: \ne.g. write \'W\' rather than \'w\', or \n\'Pb\' rather than \'pb\' " << endl;
                cerr << "\nABORT CALCULATION" << endl;
                exit(1);
            } else
                cout << "Gamma-ray absorption coefficients for " << str_pure_element[i] << " taken from file:\n"
                     << filename << endl;

            nlines = 0;

            while (true) {

                infile >> buffer.energy >> buffer.u_rho;

                in.push_back(buffer);

                if (infile.eof() || infile.fail())
                    break;

                nlines++;
            }
            in.pop_back();
            cout << "Number of lines read: " << nlines << endl;
            infile.close();

            //The next 'for' loop is needed to fill the [u_gamma_weighted] array
            double u_gamma_weighted[MAX_NUMBER_LINES] = { 0 };

            int k = 0;
            for (std::vector<record>::iterator j = in.begin(); j != in.end(); ++j) {
                //int energy(j->energy);
                double u_rho(j->u_rho);

                //energy_list[k] = energy;
                energy_list[k] = &j->energy;
                u_gamma_weighted[k] = Function.getWeightedMF(u_rho, w_mass_factor);

                k++;
            }

            //The next 'for' loop should be used to do event-by-event physics
            for (int k = 0; k < nlines; k++) {
                compound_mu_gamma[k] += u_gamma_weighted[k];
            }
        }

        double u_neutron = compound_mu_neutron;
        double d_u_neutron = err_compound_mu_neutron;

        //outfile<<"Energy\tI/I0\tdI_I0"<<endl;
        outfile << setiosflags(ios::left);
        outfile << setw(15) << "Energy"
                << "\t" << setw(20) << "u_g [cm^{-1}]"
                << "\t" << setw(20) << "I/I0"
                << "\t" << setw(20) << "dI_I0" << endl;
        outfile << endl;

        //The next 'for' loop should be used to do SUMMED event-by-event physics
        for (int k = 0; k < nlines; k++) {
            double u_gamma = compound_mu_gamma[k] * Number.getCompoundDensity();
            if (k == 164)
                cout << "500 keV " << k << " " << compound_mu_gamma[k] << endl;

            double d_u_gamma = (u_gamma / 100.0) * 5.0;

            double I_I0 = Function.calcAttenuation(Number.getThickness(), u_gamma, u_neutron, alpha_rad, cmORmm, uGammaOnly);

            double d_I_I0 = Function.getErrAttenuation(Number.getThickness(), alpha_rad, I_I0, cmORmm, u_neutron, u_gamma, d_u_neutron, d_u_gamma, uGammaOnly);

            //outfile<<*energy_list[k]<<I_I0<<d_I_I0<<endl;
            outfile << setiosflags(ios::left);
            //outfile<<setiosflags(ios::fixed);
            outfile << setiosflags(ios::scientific);
            outfile << setw(15) << *energy_list[k] << "\t" << setw(20) << setprecision(9) << compound_mu_gamma[k] << "\t" << I_I0 << "\t" << setw(20) << d_I_I0 << endl;
        }
    }

    cout << endl;
    cout << "Results of self-attenuation calculations written to file." << endl;

    outfile.close();

    return 0;
}
