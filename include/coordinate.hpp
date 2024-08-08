#ifndef COORDINATES_HEADER
#define COORDINATES_HEADER
#include <cstdint>
enum coordate_type : std::uint8_t { CARTESIAN = 0, SPHERICAL, CYLINDRICAL };

class coordinate_transformer
{
public:
    static void transform( unsigned int& num, double* data, coordate_type& from,
                           coordate_type& to );

#ifdef DEBUG
public:
#else
private:
#endif
    static void car2sph( double data[ 3 ] );
    static void sph2car( double data[ 3 ] );
    static void car2cyl( double data[ 3 ] );
    static void cyl2car( double data[ 3 ] );
    static void sph2cyl( double data[ 3 ] );
    static void cyl2sph( double data[ 3 ] );
};
#endif
