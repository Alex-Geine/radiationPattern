#include "RadiationPattern.h"

#define Pi_2 (M_PI * 2)

// Constructor
Point3::Point3(double x, double y, double z) : m_x(x), m_y(y), m_z(z)
{
}

// Get z coordinate
void Point3::GetZCoord(double r)
{
    m_z = sqrt(r * r - m_x * m_x - m_y * m_y );

    return;
}

// Get distanse between two points
double Point3::GetR(Point3 p)
{
    return sqrt((m_x - p.m_x) * (m_x - p.m_x) + (m_y - p.m_y) * (m_y - p.m_y) + (m_z - p.m_z) * (m_z - p.m_z));
}

// condtructor
Receiver::Receiver(Point3 p) : m_coord(p)
{
}

// Find amplitude from emmiter
std::complex<double> Receiver::FindAmpl(Point3 p, double l)
{    
    double r = p.GetR(m_coord);

    std::complex<double> phase(0.0, r * Pi_2 / l);

    return exp(-phase) / r;
}

// Constructor
Antenna::Antenna(std::vector<std::vector<bool>> receiversLocation, double d) : m_d(d)
{
    uint64_t str = receiversLocation.size();
    uint64_t col = receiversLocation[0].size();

    double left = - (double)col / 2. * m_d + m_d / 2.;
    double bot  = - (double)str / 2. * m_d + m_d / 2.;

    double x = 0;
    double y = 0;
    double z = 0;

    for (uint64_t i = 0; i < str; ++i)
    {
        for (uint64_t j = 0; j < col; ++j)
            if (receiversLocation[i][j])
            {
                x = left + (double)j * m_d;
                y = bot + (double)i * m_d;

                m_receivers.push_back(Receiver(Point3(x, y)));
            }
    }
}

// Get Amplitude of signal
double Antenna::GetSigAmpl(Point3 p, double l)
{
    std::complex<double> sum{0,0};
    for (uint64_t i = 0; i < m_receivers.size(); ++i)
        sum += m_receivers.at(i).FindAmpl(p, l);

    return sqrt(std::abs(sum));
}

bool check_cirle(double x, double y, double r)
{
    if (x * x + y * y > r * r)
        return false;

    return true;
}

// Function for creating radiation pattern
// [in] locations - 2d array with antennas locations (true/false)
// [in] l - wave length
// [in] k - paramener for distance between receivers (d = k * l)
// [in] r - radius of area
// [in] N - number of points to compute
// return 1d array with amplitude values
double* g_getDiagram(std::vector<std::vector<bool>> locations, double l, double k, double r, int n)
{
    double* ampl = new double[n * n];

    // Init antenna
    Antenna antenna(locations, l * k);

    double dx = 2. * r / ((double)n - 1.);

    double left = -r;  // Left coordinate of square
    double bot  = left;     // Bot coordinate of square

    double x = 0;
    double y = 0;

    for (int raw = 0; raw < n; ++raw)
    {
        x = left + raw * dx;

        for (int col = 0; col < n; ++col)
        {
            y = bot + dx * col;

            if (check_cirle(x, y, r))
            {
                Point3 p (x, y);
                p.GetZCoord(r);
                ampl[raw * n + col] = antenna.GetSigAmpl(p, l);
            }
            else
            {
                ampl[raw * n + col] = 0.;
            }
        }
    }
    return ampl;
}
