#include "RadiationPattern.h"
#include "Signal.h"
#include "Picture.h"
#include <iostream>

// log scale for real data
void logScale(std::complex<double>** data, uint32_t col, uint32_t str)
{
    for (uint64_t i = 0; i < str; ++i)
        for (uint64_t j = 0; j < col; ++j)
            data[i][j] = std::log10(std::complex<double>(1,1) + data[i][j]);
}

int main()
{
    int n = 0;
    int P = 1000;
    double l = 1;
    double k = 0.5;
    double r = 1000;

    std::cout << "Type antennas: " << std::endl;
    std::cin >> n;

    std::cout << "Type points: " << std::endl;
    std::cin >> P;

    std::cout << "Antennas: " << std::endl;
    std::vector<std::vector<bool>> loc;
    loc.resize(n);
    
    for (int i = 0; i < n; ++i)
    {
        loc[i].resize(n);
        for (int j = 0; j < n; ++j)
        {
            loc[i][j] = true;
        //    std::cout << loc[i][j] << " ";
        }
      //  std::cout << std::endl;
    }
    //std::cout << std::endl;
    
    double* ampl = g_getDiagram(loc, l, k, r, P);
    std::cout << "Diagramm is ok!" << std::endl;
    /*
    std::cout << "ampl:" << std::endl;
    for (int i = 0; i < P; ++i)
    {
        for (int j = 0; j < P; ++j)
        {
           std::cout << ampl[i * P + j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    */

    Signal sig(P, P);
  
    std::complex<double>** data = sig.GetDataArray();

    for (int i = 0; i < P; ++i)
        for (int j = 0; j < P; ++j)
            data[i][j] = std::complex<double>{ampl[i * P + j], 0};

    std::cout << "Signal!" << std::endl;
    //logScale(data, P, P);

    uint8_t* dataPic = new uint8_t[P * P];
    uint8_t* picture = new uint8_t[P * P * 3];

    sig.GetPicture(dataPic, false);
    std::cout << "Picture!" << std::endl;
    g_toGrayScaleOut(P, P, dataPic, picture);

    g_safeImage(std::string("Diagramm.png"), P, P, picture);
    
    return 0;
}