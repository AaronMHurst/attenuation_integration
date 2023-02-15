#ifndef INPUTLIB_H
#define INPUTLIB_H

class record{
public:
  int energy;
  double u_rho;

  record(int i1, double d1) :
    energy(i1), u_rho(d1) {}
  record() {}
  ~record() {}
};

class physical_data{
public:
  std::string chemical_symbol;
  int atomic_number;
  double A, rho, elemental_sigmaNABS, err_elemental_sigmaNABS;
  double incoh_scatter_cs, err_incoh_scatter_cs;

  physical_data(std::string s1, int i1, double d1, double d2, double d3, double d4, double d5, double d6) : chemical_symbol(s1), atomic_number(i1), A(d1), rho(d2), elemental_sigmaNABS(d3), err_elemental_sigmaNABS(d4), incoh_scatter_cs(d5), err_incoh_scatter_cs(d6) {}
  physical_data() {}
  ~physical_data() {}
};

#endif

