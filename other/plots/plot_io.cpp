#include <string>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

// ------------ LÊ CSV COM 2 COLUNAS, PARA N 
void read_two_columns(
    const std::string& filename,
    std::vector<double>& x,
    std::vector<double>& y)
{
    std::ifstream file(filename);
    std::cout << "Reading file: " << filename << std::endl;
    if(!file)
        throw std::runtime_error("Cannot open file: " + filename);

    std::string line;

    std::getline(file,line); // header

    while(std::getline(file,line))
    {
        if(line.empty()) continue;

        std::stringstream ss(line);
        std::string a,b;

        std::getline(ss,a,',');
        std::getline(ss,b,',');

        x.push_back(std::stod(a));
        y.push_back(std::stod(b));
    }
}

void read_sigma_exp(
    const std::string& filename,
    std::vector<int>& dataset,
    std::vector<double>& W,
    std::vector<double>& sigma,
    std::vector<double>& err)
{
    std::ifstream file(filename);
    std::string line;

    while (std::getline(file, line))
    {
        // remove espaços iniciais
        line.erase(0, line.find_first_not_of(" \t"));

        if (line.empty() || line[0] == '#')
            continue;

        std::vector<double> cols;

        // --- tenta parsing como CSV ---
        {
            std::stringstream ss(line);
            std::string token;

            while (std::getline(ss, token, ','))
            {
                try {
                    cols.push_back(std::stod(token));
                }
                catch (...) {
                    // ignora tokens não numéricos (ex: "-")
                }
            }
        }

        // --- fallback: separação por espaço/tab ---
        if (cols.size() < 3)
        {
            cols.clear();
            std::stringstream ss(line);
            std::string token;

            while (ss >> token)
            {
                try {
                    cols.push_back(std::stod(token));
                }
                catch (...) {
                    // ignora lixo
                }
            }
        }

        if (cols.size() < 3)
            continue;

        // --- caso 1: formato simples ---
        if (cols.size() == 3)
        {
            int d = 0;
            dataset.push_back(d);
            W.push_back(cols[0]);
            sigma.push_back(cols[1]);
            err.push_back(cols[2]);
        }

        // --- caso 2: erro assimétrico (4 colunas) ---
        else if (cols.size() >= 4 && cols.size() < 8)
        {
            int d = 0;

            double W_val = cols[0];
            double s_val = cols[1];
            double dy_minus = std::abs(cols[2]);
            double dy_plus  = std::abs(cols[3]);

            double e_val = 0.5 * (dy_minus + dy_plus);

            dataset.push_back(d);
            W.push_back(W_val);
            sigma.push_back(s_val);
            err.push_back(e_val);
        }

        // --- caso 3: HEPData completo ---
        else if (cols.size() >= 8)
        {
            int d = static_cast<int>(cols[0]);

            double W_val     = cols[1];
            double sigma_val = cols[3];

            double stat_p = cols[4];
            double stat_m = std::abs(cols[5]);
            double sys_p  = cols[6];
            double sys_m  = std::abs(cols[7]);

            double stat = 0.5 * (stat_p + stat_m);
            double sys  = 0.5 * (sys_p + sys_m);

            double err_val = std::sqrt(stat*stat + sys*sys);

            dataset.push_back(d);
            W.push_back(W_val);
            sigma.push_back(sigma_val * 1000.0); // μb → nb
            err.push_back(err_val * 1000.0);
        }
    }
}


//-------------- LEITOR DE DADOS EXPERIMENTAIS DO XYSCAN ----------------------
void read_xyscan_csv(const std::string& filename,
                     std::vector<double>& W,
                     std::vector<double>& sigma,
                     std::vector<double>& error)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Erro ao abrir arquivo: " << filename << "\n";
        return;
    }

    std::string line;

    while (std::getline(file, line)) {

        // ignora comentários e linhas vazias
        if (line.empty() || line[0] == '#')
            continue;

        std::stringstream ss(line);
        std::string col;
        std::vector<double> values;

        while (std::getline(ss, col, ',')) {
            if (col.empty()) continue; // ignora vírgula final
            values.push_back(std::stod(col));
        }

        if (values.size() < 3)
    continue;

double w   = values[0];
double sig = values[1];
double err;

// caso 1: CSV com 3 colunas (erro já simétrico)
if (values.size() == 3) {
    err = values[2];
}
// caso 2: CSV com 4 colunas (erro assimétrico)
else {
    double dy_minus = values[2];
    double dy_plus  = values[3];
    err = 0.5 * (dy_minus + dy_plus);
}

W.push_back(w);
sigma.push_back(sig);
error.push_back(err);
    }

    file.close();
}

// ------------ LÊ CSV'S COM 1 INDEPENDENTE E DUAS DEPENDENTES, IDEAL PARA COMPARAÇÕES
// ------------- BOOSTED GAUSSIAN - GAUSSIAN LIGHT CONE

void read_csv(
    const std::string& filename,
    std::vector<double>& x,
    std::vector<double>& glc,
    std::vector<double>& bg)
{
    std::ifstream file(filename);

    if(!file)
        throw std::runtime_error("Cannot open file: " + filename);

    x.clear();
    glc.clear();
    bg.clear();

    std::string line;

    std::getline(file, line); // header

    const double MIN_PLOT_VALUE = 1e-8;  // Threshold to avoid matplotlib numerical issues

    while (std::getline(file, line))
    {
        if(line.empty()) continue;

        std::stringstream ss(line);
        std::string a,b,c;

        if (!std::getline(ss, a, ',')) continue;
        if (!std::getline(ss, b, ',')) continue;
        if (!std::getline(ss, c, ',')) continue;

        try {
            double x_val = std::stod(a);
            double glc_val = std::stod(b);
            double bg_val = std::stod(c);

            // Filter out unreasonably small values that cause matplotlib-cpp issues
            if (std::abs(glc_val) < MIN_PLOT_VALUE) glc_val = MIN_PLOT_VALUE;
            if (std::abs(bg_val) < MIN_PLOT_VALUE) bg_val = MIN_PLOT_VALUE;

            // Check for NaN/Inf
            if (!std::isfinite(x_val) || !std::isfinite(glc_val) || !std::isfinite(bg_val)) {
                std::cerr << "Warning: Skipping line with NaN/Inf: " << line << std::endl;
                continue;
            }

            x.push_back(x_val);
            glc.push_back(glc_val);
            bg.push_back(bg_val);
        }
        catch (...) {
            // ignora linha mal formatada
            continue;
        }
    }

    // sanity check forte
    if (x.size() != glc.size() || x.size() != bg.size()) {
        throw std::runtime_error("CSV columns size mismatch in: " + filename);
    }
}

void read_rapidity_hepdata(
    const std::string& filename,
    std::vector<double>& y,
    std::vector<double>& dsdy,
    std::vector<double>& err)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file: " << filename << std::endl;
        return;
    }

    std::string line;

    while (std::getline(file, line))
    {
        if (line.empty()) continue;

        //  remove espaços iniciais
        line.erase(0, line.find_first_not_of(" \t"));

        if (line.empty()) continue;

        // ignora headers
        if (line[0] == '#' || line[0] == '|')
            continue;
        std::stringstream ss(line);
        std::string token;
        std::vector<double> cols;

        while (std::getline(ss, token, ','))
        {
            try {
                cols.push_back(std::stod(token));
            } catch (...) {}
        }

        if (cols.size() < 8) continue;

        double y_val = cols[0];
        double sigma = cols[3];

        double stat_p = cols[4];
        double stat_m = std::abs(cols[5]);
        double sys_p  = cols[6];
        double sys_m  = std::abs(cols[7]);

        double stat = 0.5 * (stat_p + stat_m);
        double sys  = 0.5 * (sys_p + sys_m);

        double error = std::sqrt(stat*stat + sys*sys);

        y.push_back(y_val);
        dsdy.push_back(sigma);
        err.push_back(error);
    }

    file.close();
    std::cout << "read_rapidity_hepdata: loaded " << y.size() << " points from " << filename << std::endl;
}