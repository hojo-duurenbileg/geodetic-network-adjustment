// Geodetic network adjustment: example driver.
//
// Reads control/approximate coordinates and horizontal angle and distance
// observations, runs a parametric least-squares adjustment, and reports the
// adjusted coordinates, variance-covariance matrix and error ellipses.
//
// Author: Khongorzul Duurenbileg, September 2019

#include <iostream>
#include <iomanip>

#include "GeoNetwork.h"


using namespace std;

int main(int argc, char* argv[])
{
	const string apriori_file = argc > 1 ? argv[1] : "data/apriori.txt";
	const string obs_file = argc > 2 ? argv[2] : "data/observations.txt";


	cout << fixed << setprecision(4);
	cout.width(10);

	GeoNetwork Test;
	Test.ReadData(apriori_file);
	Test.ReadData(obs_file);
	Test.SetApriori(pow(1.5, 2));
	Test.AugmentClnP();
	Test.AugmentAnW();
	Test.DispInput();
	Test.DispCl();
	Test.DispA();
	Test.DispW();
	Test.ParaLS();
	Test.ErrorEllipses();
	Test.DispOutput();


    return 0;
}

