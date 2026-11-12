#ifndef GEONETWORKS_H
#define GEONETWORKS_H

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Includes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#include "PositioningNetwork.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <map>
#include <list>
#include <iterator>
#include <ctime>

#include <Eigen/Dense>


using namespace std;
using namespace Eigen;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Macros ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define SHOW(a) cout << #a << ": " << (a) << endl

#define HEADING(a) cout << "		" << a << endl << "		==========================================\n\n"

class GeoNetwork
{
public:
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	GeoNetwork();

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Destructors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	~GeoNetwork();

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Structs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	struct point
	{
		char type_of_point;
		double x, y, z;
		double xo, yo, zo;
		double sx, sy, sz;
		double semi_major;
		double semi_minor;
		int beta;
		int phi_max;
		int phi_min;

		void Meter2Centimeter() {
			x *= 100;
			y *= 100;
			z *= 100;
			xo *= 100;
			yo *= 100;
			zo *= 100;
			sx *= 100;
			sy *= 100;
			sz *= 100;
		}
		void Centimeter2Meter() {
			x /= 100;
			y /= 100;
			z /= 100;
			xo /= 100;
			yo /= 100;
			zo /= 100;
			sx /= 100;
			sy /= 100;
			sz /= 100;
		}
	};	
	
	struct ha
	{
		char from, at, to;
		char type_of_obs;
		double d, m, s, ss;
		double dd, sdd;
		double rads, srads;
		void DMS2DD2Radians() {
			dd = d + m / 60 + s / 3600;
			sdd = ss / 3600;
			rads = dd * acos(-1) / 180;
			srads = sdd * acos(-1) / 180;
		}
	};
	
	struct hd
	{
		char at, to;
		char type_of_obs;
		double d;
		double sd;	// centimeters
		void Meter2Centimeter() {
			d *= 100;
			sd *= 100;
		}
		void Centimeter2Meter() {
			d /= 100;
			sd /= 100;
		}
	};
		
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ File IO's ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	void ReadData(string filename);

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Mutators ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	void ParaLS();

	void AugmentAnW();

	void AugmentClnP();

	void SetApriori(double ap);
		
	void SetPoints(string filename);

	void SetObservables(string filename);

	void ErrorEllipses();
	
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Accessors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	void DispInput();

	void DispA();

	void DispW();

	void DispCl();

	void DispOutput();

private:
	
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ File IO's ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	void ReadApriori(ifstream& infile);

	void ReadControlPoint(ifstream& infile);

	void ReadAprioriOfPoint(ifstream& infile);

	void ReadHorizontalAngle(ifstream& infile);

	void ReadHorizontalDistance(ifstream& infile);

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Mutators ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	void Count();

	void AdjustX();

	void AdjustX(MatrixXd x_hat);

	void UpdateL();

	double InverseBearing(char from, char to);

	void Radians2DMS(double angle_rad);

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Constant Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	const double pi = acos(-1);		// pi
	const double rho = 3600 * 180 / pi;

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Input Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	string poe_filename;
	string obs_filename;

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	int type_of_network;	// 1 = 1D, 2 = 2D, 3 = 3D
	int num_hd_obs;
	int num_ha_obs;
	int num_va_obs;
	int u, n;
	int num_control_points;
	int num_unknown_points;
	int iterations;
	double apriori;
	double aposteriori;

	MatrixXd A;			// design matrix
	MatrixXd w;			// misclosure vector
	MatrixXd l;			// measurements vector
	MatrixXd l_adj;			// measurements vector
	MatrixXd r;			// residuals vector
	MatrixXd Cl;		// covariance matrix of observations
	MatrixXd Cl_adj;	// adjusted covariance matrix of observations
	MatrixXd P;			// weight matrix
	MatrixXd Cx;		// covariance matrix of observables
	MatrixXd Cv;		// covariance matrix of residuals
	MatrixXd x;			// best estimate of observables
	MatrixXd dx;		// correction value
	MatrixXd dists;		// horizontal distance measurements between points, format for row -->	[from, to, distance in meters] = [1, 2, 1100.200] 
	MatrixXd thetas;	// horizontal internal angular measurements between points, format for row -->	[from, at, to, theta in radians] = [1, 5, 2, 0.575684] 


	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Containers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	map<char, point> pts;				// Point container
	map<char, point>::iterator it;		// Point iterator
	map<int, ha> ha_obs;				// Number of points container
	map<int, ha>::iterator iha;			// Number of points iterator
	map<int, hd> hd_obs;				// Number of points container
	map<int, hd>::iterator ihd;			// Number of points iterator

};


#endif // !GEONETWORKS_H