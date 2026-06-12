#include "../other/integration.hpp"
#include "nuclear.hpp"
#include "../other/ctes.hpp"

#include <fstream>
#include <cmath>
#include <filesystem>

Nucleus Pb208("Pb", 82.0,                      // Número atômico 
                     6.62 * CFAC,               // Raio nuclear fm -> GeV^-1
                     0.546 * CFAC,              // diffusividade fm -> GeV^-1
                     0.1603/std::pow(CFAC, 3)); // densidade central fm^-3 -> GeV^3
                                                // parâmetros típicos para Pb-208

Nucleus proton("p", 1.0,                      // Número atômico
                     0.84 * CFAC,              // Raio do próton fm -> GeV^-1
                     0.0,                      // sem diffusividade
                     0.17/std::pow(CFAC, 3)); // densidade central ajustada para A=1 fm^-3 -> GeV^3
                                                // parâmetros típicos para próton

TA_Table precompute_TA(int Nb, double bmax, const Nucleus& N)
{
    TA_Table table;
    table.b_vals.resize(Nb);
    table.TA_vals.resize(Nb);

    double zmax = N.R * 5.07; // fm para GeV^-1, suficiente para Pb-208
    int Nz = 200;

    for (int i = 0; i < Nb; ++i)
    {
        double b = bmax * i / (Nb - 1);
        table.b_vals[i] = b;


        std::function<double(double)> integrand = [b, &N](double z) {
            double r = sqrt(b*b + z*z);
            return rho_WS(r, N);
        };

        table.TA_vals[i] = integrate_simpson(integrand, -zmax, zmax, Nz);
    }

    return table;
}



double interpolate_TA(double b, const TA_Table& table)
{
    int N = table.b_vals.size();

    if (b <= table.b_vals[0]) return table.TA_vals[0];
    if (b >= table.b_vals[N-1]) return table.TA_vals[N-1];
    // busca binária
    int low = 0, high = N - 1;
    while (high - low > 1)
    {
        int mid = (low + high) / 2;
        if (table.b_vals[mid] > b)
            high = mid;
        else
            low = mid;
    }

    double b1 = table.b_vals[low];
    double b2 = table.b_vals[high];
    double T1 = table.TA_vals[low];
    double T2 = table.TA_vals[high];

    return T1 + (T2 - T1) * (b - b1) / (b2 - b1);
}

double rho_WS(double r, const Nucleus& N)
{
    return N.rho0 / (1.0 + exp((r - N.R)/N.a));
}

double integral_rho(const Nucleus& N)
{
    double rmax = 20.0; // fm (suficiente para Pb)
    int Nr = 500;

    auto integrand = [&N](double r) {
        double rho = 1.0 / (1.0 + exp((r - N.R)/N.a)); // rho0 = 1
        return r*r * rho;
    };

    double integral = integrate_simpson(integrand, 0.0, rmax, Nr);

    return 4.0 * M_PI * integral;
}

double compute_rho0(int A, const Nucleus& N)
{
    double I = integral_rho(N);
    return A / I;
}

namespace fs = std::filesystem;


// --------------------------------------------------
// salva tabela
// --------------------------------------------------
void save_TA_table(const TA_Table& table,
                   const std::string& filename)
{
    std::ofstream fout(filename);

    if (!fout)
    {
        throw std::runtime_error(
            "Erro ao abrir arquivo para escrita: "
            + filename);
    }

    fout << "# b[GeV^-1] TA(b)\n";

    for (size_t i = 0; i < table.b_vals.size(); ++i)
    {
        fout << table.b_vals[i]
             << " "
             << table.TA_vals[i]
             << "\n";
    }

    fout.close();
}

// --------------------------------------------------
// carrega tabela
// --------------------------------------------------
TA_Table load_TA_table(const std::string& filename)
{
    TA_Table table;

    std::ifstream fin(filename);

    if (!fin)
    {
        throw std::runtime_error(
            "Erro ao abrir arquivo: "
            + filename);
    }

    std::string line;

    while (std::getline(fin, line))
    {
        if (line.empty() || line[0] == '#')
            continue;

        std::stringstream ss(line);

        double b, TA;

        ss >> b >> TA;

        table.b_vals.push_back(b);
        table.TA_vals.push_back(TA);
    }

    fin.close();

    return table;
}

// --------------------------------------------------
// função principal
// --------------------------------------------------
TA_Table get_TA_table(const Nucleus& nuc,
                      int Nb,
                      double bmax)
{
    // cria diretório se não existir
    fs::create_directories("out/csv/Ta_tables");

    std::string filename =
        "out/csv/Ta_tables/TA_" + nuc.name + ".dat";

    // -----------------------------------------
    // se arquivo existir -> carrega
    // -----------------------------------------
    if (fs::exists(filename))
    {
        std::cout << "Carregando TA_table de "
                  << filename << "\n";

        return load_TA_table(filename);
    }

    // -----------------------------------------
    // senão -> calcula
    // -----------------------------------------
    std::cout << "TA_table não encontrada para "
              << nuc.name
              << ". Gerando...\n";

    TA_Table table =
        precompute_TA(Nb, bmax, nuc);

    // salva
    save_TA_table(table, filename);

    std::cout << "TA_table salva em "
              << filename
              << "\n";

    return table;
}