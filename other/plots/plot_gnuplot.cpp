#include <gnuplot-iostream.h>
#include <vector>
#include <utility>

void teste(){
    Gnuplot gp;

    std::vector<std::pair<double,double>> pts;

    for(int i=0;i<100;i++) {
        double x = i*0.1;
        double y = sin(x);

        pts.emplace_back(x,y);
    }

    gp << "plot '-' with lines title 'sin(x)'\n";

    gp.send1d(pts);
}