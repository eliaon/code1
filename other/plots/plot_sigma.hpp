#ifndef PLOT_SIGMA_HPP
#define PLOT_SIGMA_HPP

#include <iostream>
#include "../../calculations/nuclear.hpp"

void plot_dsigma_dt(std::string meson);

void plot_sigma(std::string meson, std::string csv_file, std::string dipolemodel = "GBW");

void plot_sigma_phi(std::string csv_file, std::string dipolemodel);
void plot_sigma_Jpsi(std::string csv_file, std::string dipolemodel);

void plot_sigma_Jpsi_gammaPb(std::string csv_file, std::string dipolemodel);
void plot_sigma_phi_gammaPb(std::string csv_file, std::string dipolemodel);
void plot_sigma_placeholder(std::string csv_file, std::string dipolemodel, std::string meson);
void plot_sigma_placeholder(std::string csv_file, std::string dipolemodel, std::string meson, Nucleus& nucleus);

void plot_sigma_gammaA(std::string meson, std::string csv, std::string modelo, Nucleus& nucleus);
void plot_sigma_models(const std::vector<std::pair<std::string,std::string>>& files, const Meson& M,
                   bool fc);











    
#endif
