#ifndef PLOT_OTHER_HPP
#define PLOT_OTHER_HPP

#include <iostream>
#include <vector>
#include "../utils.hpp"

void plot_TA_b();

void compare_N_models(double x);

void plot_overlap(std::string meson);

void plot_N_dglap(std::string csv_file1, std::string csv_file2);

void plot_overlap_fc(std::string csv_file, std::string csv_file_fc, std::string meson);

void plot_xf_multiQ2();
void plot_N_multi(
    const std::vector<std::string>& filenames,
    const std::vector<double>& x_vals,
    const std::string& modelo
);

void plot_N_models(const std::vector<std::pair<std::string,std::string>>& files, double x);

void plot_XY(const std::vector<double>& x,
             const std::vector<double>& y,
             const std::string& xlabel,
             const std::string& ylabel,
             const std::string& title);
             
void plot_overlap_models(std::vector<std::pair<std::string, std::string>>& files, const Meson& M, bool fc);


#endif