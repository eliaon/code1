#ifndef PLOT_RAPIDITY_HPP
#define PLOT_RAPIDITY_HPP

#include <iostream>


void plot_rapidity();

void plot_rapidez_PbPb_Jpsi(std::string csv_file, std::string model, double sqrt_s);
void plot_rapidez_PbPb_phi(std::string csv, std::string model, double sqrt_s);
void plot_rapidez(std::string csv, std::string model, std::string meson, std::string collision,double sqrt_s);


#endif