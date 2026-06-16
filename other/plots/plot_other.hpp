#ifndef PLOT_OTHER_HPP
#define PLOT_OTHER_HPP

#include <iostream>
#include <vector>

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

void plot_XY(std::vector<double> x, std::vector<double> y, std::string xlabel, std::string ylabel, std::string title);



#endif