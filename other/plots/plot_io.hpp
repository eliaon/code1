#ifndef PLOT_IO_HPP
#define PLOT_IO_HPP


#include <iostream>
#include <vector>
#include <string>


void read_two_columns(
    const std::string& filename,
    std::vector<double>& x,
    std::vector<double>& y);

void read_sigma_exp(
    const std::string& filename,
    std::vector<int>& dataset,
    std::vector<double>& W,
    std::vector<double>& sigma,
    std::vector<double>& err);

void read_xyscan_csv(const std::string& filename,
                     std::vector<double>& W,
                     std::vector<double>& sigma,
                     std::vector<double>& error);

void read_csv(
    const std::string& filename,
    std::vector<double>& x,
    std::vector<double>& glc,
    std::vector<double>& bg);

void read_rapidity_hepdata(
    const std::string& filename,
    std::vector<double>& y,
    std::vector<double>& dsdy,
    std::vector<double>& err);


#endif // PLOT_IO_HPP