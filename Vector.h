//https://apps.dtic.mil/sti/pdfs/AD1023343.pdf
// Vector.h: Definition & implementation of class for the algebra of 3D vectors
// R. Saucier, February 2000 (Revised June 2016)
#ifndef VECTOR_3D_H
#define VECTOR_3D_H
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <iostream>

namespace va { // vector algeba namespace
enum rep { CART, POLAR }; // vector representation (cartesian or polar)
enum comp { X, Y, Z }; // for referencing cartesian components
const double R2D( 180. / M_PI ); // for converting radians to degrees
const double D2R( M_PI / 180. ); // for converting degrees to radians
inline double deg( const double rad ) { return rad * R2D; } // convert radians to degrees
inline double rad( const double deg ) { return deg * D2R; } // convert degrees to radians
// IEEE 754 standard has 53 significant bits or 15 decimal digits of accuracy, so anything smaller is not significant
static const long double TOL = 1.e-15;
class Vector {
// friends list
// overloaded arithmetic operators

	friend Vector operator+( const Vector& a, const Vector& b ) { // addition
 		return Vector( a._x + b._x, a._y + b._y, a._z + b._z );
	}

 	friend Vector operator-( const Vector& a, const Vector& b ) { // subtraction
 		return Vector( a._x - b._x, a._y - b._y, a._z - b._z );
	}

	friend Vector operator*( const Vector& v, double s ) { // right multiply by scalar
		return Vector( s * v._x, s * v._y, s * v._z );
	}

	friend Vector operator*( double s, const Vector& v ) { // left multiply by scalar
		return Vector( s * v._x, s * v._y, s * v._z );
	}

	friend Vector operator/( const Vector& v, double s ) { // division by scalar
		assert( s != 0. );
		return Vector( v._x / s, v._y / s, v._z / s );
	}
	
	// overloaded algebraic operators
	friend inline double operator*( const Vector& a, const Vector& b ) { // dot product
		double c( a._x * b._x +
				  a._y * b._y +
				  a._z * b._z );
		if ( fabs(c) < TOL ) c = 0.0L; // set precision to be no more than 15 digits
		return c;
 	}

 	friend inline Vector operator^( const Vector& a, const Vector& b ) { // cross product
		return Vector( a._y * b._z - a._z * b._y,
					   a._z * b._x - a._x * b._z,
					   a._x * b._y - a._y * b._x );
	}

 	// access functions
	friend double x( const Vector& v ) { // x-coordinate
 		return v._x;
	}

	friend double y( const Vector& v ) { // y-coordinate
 		return v._y;
	}
	
	friend double z( const Vector& v ) { // z-coordinate
		return v._z;
	}

	friend double r( const Vector& v ) { // magnitude of vector
		return v._mag();
	}

	friend double theta( const Vector& v ) { // polar angle (radians)
		return v.theta();
	}

	friend double phi( const Vector& v ) { // azimuthal angle (radians)
		return v.phi();
	}

	friend double norm( const Vector& v ) { // norm or magnitude
		return v._mag();
	}

	friend double mag( const Vector& v ) { // magnitude
		return v._mag();
	}

	friend double scalar( const Vector& v ) { // magnitude
		return v._mag();
	}

	friend double angle( const Vector& a, const Vector& b ) { // angle (radians) between vectors

 		double s = a.unit() * b.unit();
 		if ( s >= 1. )
 			return 0.;
 		else if ( s <= -1. )
 			return M_PI;
 		else
 		return acos( s );
 	}

 	friend Vector unit( const Vector& v ) { // unit vector in same direction
		return v.unit();
	}

	friend Vector normalize( const Vector& v ) { // returns a unit vector in same direction
		return v.unit();
	}

 	friend double dircos( const Vector& v, const comp& i ) { // direction cosine
		return v.dircos( i );
	}
	friend Vector proj( const Vector& a, // projection of first vector
	 					const Vector& b ) { // along second vector
		return a.proj( b );
	}


	// overloaded stream operators
	friend std::istream& operator>>( std::istream& is, Vector& v ) { // input vector
		return is >> v._x >> v._y >> v._z;
	}

	friend std::ostream& operator<<( std::ostream& os, const Vector& v ) { // output vector
		Vector a( v );
		a._set_precision();
		return os << a._x << " " << a._y << " " << a._z;
	}

public:
	Vector( double x, double y, double z, // constructor (cartesian or polar
			rep mode = CART ) { // with cartesian as default)
 		if ( mode == CART ) { // cartesian form
			this->_x = x;
			this->_y = y;
			this->_z = z;
		}
		else if ( mode == POLAR ) // polar form
			_setCartesian( x, y, z );
		else {
			std::cerr << "Vector: mode must be either CART or POLAR" << std::endl;
			exit( EXIT_FAILURE );
		}
	}

	// Vector( std::array<double,3> a, rep mode): Vector(a[0], a[1], a[2], mode){};

	Vector( void ) : _x( 0. ), _y( 0. ), _z( 0. ) { // default constructor
	}

	~Vector( void ) { // default destructor
	}

	Vector( const Vector& v ) : _x( v._x ), // copy constructor
								_y( v._y ),
								_z( v._z ) {
 	}

	Vector& operator=( const Vector& v ) { // assignment operator
		if ( this != &v ) {
			_x = v._x;
			_y = v._y;
			_z = v._z;
		}
		return *this;
	}

	// overloaded arithmetic operators
	Vector& operator+=( const Vector& v ) { // addition assignment
 		_x += v._x;
 		_y += v._y;
 		_z += v._z;
 		return *this;
 	}
	Vector& operator-=( const Vector& v ) { // subtraction assignment
		_x -= v._x;
		_y -= v._y;
		_z -= v._z;
		return *this;
	}

	Vector& operator*=( double s ) { // multiplication assignment

		_x *= s;
		_y *= s;
		_z *= s;
		return *this;
	}

	Vector& operator/=( double s ) { // division assignment

		assert( s != 0. );
		_x /= s;
		_y /= s;
		_z /= s;
		return *this;
	}

	Vector operator-( void ) { // negative of a vector

		return Vector( -_x, -_y, -_z );
	}

	const double& operator[]( comp i ) const { // index operator (component)
		if ( i == X ) return _x;
		if ( i == Y ) return _y;
		if ( i == Z ) return _z;
		std::cerr << "Vector: Array index out of range; must be X, Y, or Z" << std::endl;
		exit( EXIT_FAILURE );
	}

	// access functions
	double x( void ) const { // x-component
		return _x;
	}

	double y( void ) const { // y-component
		return _y;
	}

	double z( void ) const { // z-component
 		return _z;
	}

	double r( void ) const { // magnitude
		return _mag();
	}

	double theta( void ) const { // polar angle (radians)
		return _theta();
	}

	double phi( void ) const { // azimuthal angle (radians)
 		return _phi();
	}

	double norm( void ) const { // norm or magnitude
		return _mag();
	}

	double mag( void ) const { // magnitude
		return _mag();
	}
 
 	operator double( void ) const { // conversion operator to return magnitude
 		return _mag();
 	}

	// utility functions
	Vector unit( void ) const { // returns a unit vector
		double m = _mag();
		if ( m > 0. ) return Vector( _x / m, _y / m, _z / m );
		std::cerr << "Vector: Cannot make a unit vector from a null vector" << std::endl;
 		exit( EXIT_FAILURE );
 	}

	Vector normalize( void ) const { // synonym for unit
		return unit();
	}

	double dircos( const comp& i ) const { // direction cosine
		if ( i == X ) return _x / _mag();
		if ( i == Y ) return _y / _mag();
		if ( i == Z ) return _z / _mag();
		std::cerr << "Vector dircos: comp out of range; must be X, Y, or Z" << std::endl;
		exit( EXIT_FAILURE );
	}

	Vector proj( const Vector& e ) const { // projects onto the given vector
		Vector u = e.unit();
		return u * ( *this * u );
	}

	Vector& rotate( const Vector& a, double angle ) { // rotates vector about given axial vector,givenangle
		Vector a_hat = a.unit(); // unit vector along the axial vector
		Vector cross = a_hat ^ *this; // cross product of a_hat with the given vector
		return *this += sin( angle ) * cross + ( 1. - cos( angle ) ) * ( a_hat ^ cross );
	}

 	Vector& rot( const Vector& a, double angle ) { // synonym for rotate
 		return rotate( a, angle );
 	}

private:

	double _x, _y, _z; // cartesian representation

	double _mag( void ) const { // compute the magnitude
		double mag = sqrt( _x *_x + _y *_y + _z *_z );
 		if ( mag < TOL ) mag = 0.0L;
 		return mag;
 	}

 	double _theta( void ) const { // compute the polar angle (radians)
		if ( _z != 0. ) return atan2( sqrt( _x * _x + _y * _y ), _z );
	 	else return M_PI_2;
	}

	double _phi( void ) const { // compute the azimuthal angle (radians)
		if ( _x != 0. )
			return atan2( _y, _x );
		else {
			if ( _y > 0 )
				return M_PI_2;
			else if ( _y == 0. )
				return 0.;
			else
				return -M_PI_2;
		}
 	}

 	void _setCartesian( double r, double theta, double phi ) {
		_x = r * sin( theta ) * cos( phi );
		_y = r * sin( theta ) * sin( phi );
		_z = r * cos( theta );
	}

	Vector _set_precision( void ) { // no more than 15 digits of precision
		if ( fabs(_x) < TOL ) _x = 0.0L;
		if ( fabs(_y) < TOL ) _y = 0.0L;
		if ( fabs(_z) < TOL ) _z = 0.0L;
 		return *this;
	}
};
	// declaration of friends
	Vector operator+( const Vector& a, const Vector& b );
	Vector operator-( const Vector& a, const Vector& b );
	Vector operator*( const Vector& v, double s );
	Vector operator*( double s, const Vector& v );
	Vector operator/( const Vector& v, double s );
	double operator*( const Vector& a, const Vector& b );
	Vector operator^( const Vector& a, const Vector& b );
	double x( const Vector& v );
	double y( const Vector& v );
	double z( const Vector& v );
	double r( const Vector& v );
	double theta( const Vector& v );
	double phi( const Vector& v );
	double mag( const Vector& v );
	double scalar( const Vector& v );
	double angle( const Vector& a, const Vector& b );
	Vector unit( const Vector& v );
	Vector normalize( const Vector& v );
	double dircos( const Vector& v, const comp& i );
	Vector proj( const Vector& a, const Vector& b );
	std::istream& operator>>( std::istream& is, Vector& v );
	std::ostream& operator<<( std::ostream& os, const Vector& v );
} // vector algebra namespace
#endif