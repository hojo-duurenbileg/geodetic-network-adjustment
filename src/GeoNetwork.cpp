#include "GeoNetwork.h"

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GeoNetwork::GeoNetwork()
{
	cout << "\n\nStarting Geodetic Network Program!!!\n\n";

	time_t rawtime;

	time(&rawtime);

	printf("Program started at local time: %s\n\n", ctime(&rawtime));

	num_hd_obs = 0;
	num_ha_obs = 0;
	num_va_obs = 0;
	num_control_points = 0;
	num_unknown_points = 0;
	iterations = 0;
	aposteriori = 1;
	apriori = 1;

}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Destructors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GeoNetwork::~GeoNetwork()
{
	cout << "Ending Geodetic Network Program\n\n";
	
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ File IO's ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GeoNetwork::ReadData(string filename) {

	cout << "Reading file: " << filename << endl << endl;

	ifstream infile;
	infile.open(filename, fstream::in);
	if (!infile.good()) {
		cout << "	File " << filename << " was not good. Try again...\n\n";
		infile.clear();  cin.ignore(256, '\n'); cin.get(); return;
	}
	char obs_type;

	while (!infile.eof()) {

		if (infile.peek() == '!') {
			poe_filename = filename;
			ReadApriori(infile);
		}
		else {
			obs_filename = filename;
			infile >> obs_type;
			switch (obs_type)
			{
			case 'H': ReadHorizontalAngle(infile); break;
			case 'D': ReadHorizontalDistance(infile); break;
			default: if (infile.eof()) { break; }
			}
		}
	}

	cout << "	Read Success: file " << filename << endl << endl;

	infile.close();
}

void GeoNetwork::ReadApriori(ifstream& infile) {
	
	infile.get();							// remove !
	infile >> type_of_network;				// type of network
	char type_of_observable;
	while (!infile.eof()) {
		if (infile >> type_of_observable) {

			switch (type_of_observable)
			{
			case 'C': ReadControlPoint(infile); break;
			case 'A': ReadAprioriOfPoint(infile); break;
			default: if (infile.eof()) { break; }
			}
		}
	}
	
}

void GeoNetwork::ReadControlPoint(ifstream& infile) {

	num_control_points++;

	char hold_id;
	point hold_point;
	hold_point.type_of_point = 'C';
	if (type_of_network == 2) 
			infile >> hold_id >> hold_point.x >> hold_point.y;
		
	hold_point.xo = hold_point.x;
	hold_point.yo = hold_point.y;
	hold_point.Meter2Centimeter();
	pts[hold_id] = hold_point;
}

void GeoNetwork::ReadAprioriOfPoint(ifstream& infile) {

	num_unknown_points++;

	char hold_id;
	point hold_point;
	hold_point.type_of_point = 'A';
	if (type_of_network == 2)
		infile >> hold_id >> hold_point.xo >> hold_point.yo;
		
	hold_point.x = hold_point.xo;
	hold_point.y = hold_point.yo;
	hold_point.Meter2Centimeter();
	pts[hold_id] = hold_point;

}

void GeoNetwork::ReadHorizontalAngle(ifstream& infile) {


	num_ha_obs++;

	ha hold_ha;
	hold_ha.type_of_obs = 'H';
	infile >> hold_ha.from >> hold_ha.at >> hold_ha.to >> hold_ha.d >> hold_ha.m >> hold_ha.s >> hold_ha.ss;

	hold_ha.DMS2DD2Radians();

	ha_obs[num_ha_obs] = hold_ha;

}

void GeoNetwork::ReadHorizontalDistance(ifstream& infile) {
	
	num_hd_obs++;

	hd hold_hd;
	hold_hd.type_of_obs = 'D';
	infile >> hold_hd.at >> hold_hd.to >> hold_hd.d >> hold_hd.sd;
	hold_hd.Meter2Centimeter();
	hd_obs[num_hd_obs] = hold_hd;

}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Mutators ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GeoNetwork::ParaLS() {

	AdjustX();
	UpdateL();

	dx = -(A.transpose() * P * A).inverse() * (A.transpose() * P * w);
	
	AdjustX(x);

	iterations++;
	while (iterations < 50 && abs(dx.maxCoeff()) > 1e-8) {

		AugmentAnW();

		dx = -(A.transpose() * P * A).inverse() * A.transpose() * P * w;

		x = x + dx;

		AdjustX(x);

		iterations++;				
	}

	r = A * dx + w;

	aposteriori = (r.transpose() * P * r).value() / (n - u);

	Cx = (A.transpose() * P * A).inverse() * aposteriori;

	l_adj = l + r;
	
	Cl_adj = A * Cx * A.transpose();
	
	Cv = aposteriori*P.inverse() - Cl_adj;
}

void GeoNetwork::AugmentAnW() {
	
	Count();
	A.setZero(n, u);
	w.setZero(n, 1);
	int i = 0, j = 0;	
	
	if (hd_obs.size() > 0 && type_of_network == 2) {

		for (ihd = hd_obs.begin(); ihd != hd_obs.end(); ihd++, i++) {
						
			char from = ihd->second.hd::at;
			char to = ihd->second.to;
			double d_from_to = sqrt(pow(pts[to].x - pts[from].x, 2) + pow(pts[to].y - pts[from].y, 2));
			
			for (j = 0, it = pts.begin(); it != pts.end(); it++) {

				if (it->second.type_of_point == 'A') {
					
					if (from == it->first) {
						A(i, j) = (pts[to].x - pts[from].x) / d_from_to;
						A(i, j + 1) = (pts[to].y - pts[from].y) / d_from_to;
					}
					else if (to == it->first) {
						A(i, j) = -(pts[to].x - pts[from].x) / d_from_to;
						A(i, j + 1) = -(pts[to].y - pts[from].y) / d_from_to;
					}

					j += 2;
				}
			}
			w(i, 0) = (ihd->second.d - d_from_to);
		}
	}

	if (ha_obs.size() > 0 && type_of_network == 2) {

		for (i = hd_obs.size(), iha = ha_obs.begin(); iha != ha_obs.end(); iha++, i++) {

			char from = iha->second.from;
			char at = iha->second.ha::at;			
			char to = iha->second.to;
			
			double d_at_from2 = pow(pts[from].x - pts[at].x, 2) + pow(pts[from].y - pts[at].y, 2);
			double d_at_to2 = pow(pts[to].x - pts[at].x, 2) + pow(pts[to].y - pts[at].y, 2);
			
			for (j = 0, it = pts.begin(); it != pts.end(); it++) {
				
				if (it->second.type_of_point == 'A') {

					if (at == it->first) {
						A(i, j) = (-(pts[from].y - pts[at].y) / d_at_from2 + (pts[to].y - pts[at].y) / d_at_to2) * rho;
						A(i, j + 1) = ((pts[from].x - pts[at].x) / d_at_from2 - (pts[to].x - pts[at].x) / d_at_to2) * rho;
					}
					else if (from == it->first) {
						A(i, j) = ((pts[from].y - pts[at].y) / d_at_from2) * rho;
						A(i, j + 1) = (-(pts[from].x - pts[at].x) / d_at_from2) * rho;
					}
					else if (to == it->first) {
						A(i, j) = (-(pts[to].y - pts[at].y) / d_at_to2) * rho;
						A(i, j + 1) = ((pts[to].x - pts[at].x) / d_at_to2) * rho;
					}

					j += 2;
				}
			}
			
			double at2from = InverseBearing(at, from);
			double at2to = InverseBearing(at, to);
			double theta_o = at2from - at2to;
			
			if (theta_o < 0.0)
				theta_o += (2 * pi);

			w(i, 0) = -(iha->second.rads - (theta_o)) * rho;

		}
	}
}

void GeoNetwork::AugmentClnP() {
	Count();
	Cl.resize(n, n);
	Cl.setZero();

	int i = 0;

	if (hd_obs.size() > 0 && type_of_network == 2) {
		for (ihd = hd_obs.begin(); ihd != hd_obs.end(); ihd++) {
			Cl(i, i) = pow(ihd->second.sd, 2);
			i++;
		}
	}
	if (ha_obs.size() > 0 && type_of_network == 2) {
		for (iha = ha_obs.begin(); iha != ha_obs.end(); iha++) {
			Cl(i, i) = pow(iha->second.ss, 2);
			i++;
		}
	}

	P = apriori * Cl.inverse();
}

void GeoNetwork::SetApriori(double ap) {
	apriori = ap;
}

void GeoNetwork::AdjustX() {

	int i = 0;
	x.resize(u, 1);
	dx.resize(u, 1);
	for (it = pts.begin(); it != pts.end(); it++) {
		if (it->second.type_of_point == 'A') {

			x(i, 0) = it->second.x;
			x(i + 1, 0) = it->second.y;
			i+=2;
		}
	}
}

void GeoNetwork::AdjustX(MatrixXd x_hat) {

	int i = 0;

	for (it = pts.begin(); it != pts.end(); it++) {
		if (it->second.type_of_point == 'A') {

			it->second.x = x_hat(i, 0);
			it->second.y = x_hat(i + 1, 0);
			i += 2;
		}
	}
}

void GeoNetwork::UpdateL() {

	r.resize(n, 1);
	l.resize(n, 1);
	int i = 0;
	
	if (hd_obs.size() > 0 && type_of_network == 2) {

		for (ihd = hd_obs.begin(); ihd != hd_obs.end(); ihd++) {
			l(i, 0) = ihd->second.d;
			i++;
		}
	}

	if (ha_obs.size() > 0 && type_of_network == 2) {

		for (iha = ha_obs.begin(); iha != ha_obs.end(); iha++) {
			l(i, 0) = iha->second.rads;
			i++;
		}
	}
}

void GeoNetwork::Count() {

	n = hd_obs.size() + ha_obs.size();
	u = num_unknown_points * type_of_network;

}

double GeoNetwork::InverseBearing(char from, char to) {

	double dx = (pts[to].x - pts[from].x);
	double dy = (pts[to].y - pts[from].y);

	// Orientate to correct quadrant
	double j = 1;
	if (dy > 0)
		j = 0;

	double t = atan(dx / dy) + pi*j;
	
	if (t < 0)
		t += (2 * pi);

	return t;
}

void GeoNetwork::Radians2DMS(double angle_rad) {

	double dd = 180 * angle_rad / pi;
	if (dd < 0) {
		cout << "-";
		dd *= -1;
	}
	else
		cout << " ";
	cout << int(dd) << "-" << int((dd - int(dd)) * 60) << "-" << setprecision(2) << ((dd - int(dd)) * 60 - int((dd - int(dd)) * 60)) * 60;
}

void GeoNetwork::SetPoints(string filename) {

}

void GeoNetwork::SetObservables(string filename) {

}

void GeoNetwork::ErrorEllipses() {

	int i = 0;
	for (it = pts.begin(); it != pts.end(); it++) {

		if (it->second.type_of_point == 'A') {

			it->second.semi_major = sqrt( 0.5 * (Cx(i, i) + Cx(i + 1, i + 1) + sqrt(pow(Cx(i, i) - Cx(i + 1, i + 1), 2) + 4 * pow(Cx(i, i + 1), 2))));

			it->second.semi_minor = sqrt(0.5 * (Cx(i, i) + Cx(i + 1, i + 1) - sqrt(pow(Cx(i, i) - Cx(i + 1, i + 1), 2) + 4 * pow(Cx(i, i + 1), 2))));


			// Orientate to correct quadrant beta
			double j = 1;
			if ((Cx(i, i) + Cx(i + 1, i + 1)) > 0)
				j = 0;

			double t = 0.5 * atan(2 * Cx(i, i + 1) / (Cx(i, i) + Cx(i + 1, i + 1))) + pi*j;
			
			if (t < 0)
				t += (2 * pi);

			it->second.beta = (int)roundl(t * 180 / pi);
		
			// Orientate to correct quadrant phi_max
			
			j = 1;
			if (Cx(i, i + 1) > 0)
				j = 0;

			t = atan(pow(it->second.semi_major, 2) - Cx(i, i) / Cx(i, i + 1));

			if (t < 0)
				t += (2 * pi);

			it->second.phi_max = (int)roundl(t * 180 / pi);

			int p2 = (int)roundl(t * 180 / pi) + 90;
			if (p2 > 360)
				p2 -= 360;
			it->second.phi_min = p2;


			i += 2;

		}
	}



}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Accessors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GeoNetwork::DispInput() {

	HEADING("Summary of Unadjusted Input Observations");
	
	if(pts.size() > 0 && type_of_network == 2) {
		HEADING("Initial a priori point estimates and fixed points (m)");
		cout << "Point	X		Y		Type\n\n";

		cout << fixed << setprecision(4);

		for (it = pts.begin(); it != pts.end(); it++) {
			it->second.Centimeter2Meter();
			cout << it->first << "	" << it->second.xo << "	" << it->second.yo << "	" << it->second.type_of_point << endl;
			it->second.Meter2Centimeter();

		}
		cout << endl;

	}

	if (hd_obs.size() > 0 && type_of_network == 2) {
		HEADING("Unadjusted Distance Observations (m)");
		cout << "From	To	Distance	StdErr	Type\n\n";

		cout << fixed << setprecision(4);

		for (ihd = hd_obs.begin(); ihd != hd_obs.end(); ihd++) {
			ihd->second.Centimeter2Meter();
			cout << ihd->second.hd::at << "	" << ihd->second.to << "	" << ihd->second.d << "	" << ihd->second.sd << "	" << ihd->second.type_of_obs << endl;
			ihd->second.Meter2Centimeter();

		}
		cout << endl;

	}

	if (ha_obs.size() > 0 && type_of_network == 2) {
		HEADING("Unadjusted Horizontal Angle Observations (DMS)");
		cout << "From	At	To	Angle		StdErr	Type\n\n";

		for (iha = ha_obs.begin(); iha != ha_obs.end(); iha++) {
			cout << fixed << setprecision(0);
			cout << iha->second.from << "	" << iha->second.ha::at << "	" << iha->second.to << "	" << iha->second.d << "-" << iha->second.m << setprecision(2) << "-" << iha->second.s << "	" << iha->second.ss << "	" << iha->second.type_of_obs << endl;
		}
		cout << endl;

	}

}

void GeoNetwork::DispA() {

	HEADING("Design Matrix");

	IOFormat CleanFmt(4, 0, " ", "\n", "", "");
	cout << "A  (row 1-6 [units 1/cm]) (row 7-18 [units \"/cm])	" << "(" << n << " x " << u << "):\n\n";
	cout << A.format(CleanFmt) << endl << endl;
}

void GeoNetwork::DispW() {

	HEADING("Misclosure Vector");

	IOFormat CleanFmt2(3, 0, " ", "\n", "", " cm");
	cout << "w " << "(" << n << " x 1" << "):\n\n";
	cout << w.block(0, 0, num_hd_obs, 1).format(CleanFmt2) << endl;

	for (int i = num_hd_obs; i < (num_ha_obs + num_hd_obs); i++) {

		cout << std::right << setw(6) << setprecision(2) << w(i, 0);
		cout << " \"" << endl;
	}
}

void GeoNetwork::DispCl() {

	HEADING("Covariance of Observations");

	IOFormat CleanFmt(2, 0, " ", "\n", "", "");
	cout << "\nCl (row-col 1-6 [units cm^2]) (row-col 7-18 [units (\")^2])	" << "(" << n << " x " << n << "):\n\n";
	cout << Cl.format(CleanFmt) << endl << endl;
}

void GeoNetwork::DispOutput() {
	
	HEADING("Summary of Adjusted Observations");
	
	if (hd_obs.size() > 0 && type_of_network == 2) {
		HEADING("Adjusted Distance Observations (meters)");
		cout << "From	To	Distance	StdErr	Type\n\n";

		cout << fixed << setprecision(4);
		int i = 0;
		for (ihd = hd_obs.begin(); ihd != hd_obs.end(); ihd++) {
			cout << ihd->second.hd::at << "	" << ihd->second.to << "	" << l_adj(i, 0)/100 << "	" << sqrt(Cl_adj(i,i)) << "	" << ihd->second.type_of_obs << endl;
			i++;
		}
		cout << endl;

	}

	if (ha_obs.size() > 0 && type_of_network == 2) {
		HEADING("Adjusted Horizontal Angle Observations (DMS)");
		cout << "From	At	To	Angle		StdErr	Type\n\n";
		int i = num_hd_obs;
		for (iha = ha_obs.begin(); iha != ha_obs.end(); iha++) {
			cout << fixed << setprecision(0);
			cout << iha->second.from << "	" << iha->second.ha::at << "	" << iha->second.to << "	";
			Radians2DMS(l_adj(i, 0));
			cout << "	" << sqrt(Cl_adj(i,i)) << "	" << iha->second.type_of_obs << endl;
			i++;
		}
		cout << endl;

	}

	HEADING("Covariance of Adjusted Observations");

	IOFormat CleanFmt4(2, 0, " ", "\n", "", "");
	cout << "Cl_adj (row-col 1-6 [units cm^2]) (row-col 7-18 [units (\")^2])	" << "(" << n << " x " << n << "):\n\n";
	cout << Cl_adj.format(CleanFmt4) << endl << endl;

	if (pts.size() > 0 && type_of_network == 2) {
		HEADING("Final point estimates and fixed points (meters)");
		cout << "Point	X		Y		Type\n\n";

		cout << fixed << setprecision(4);

		for (it = pts.begin(); it != pts.end(); it++) {
			it->second.Centimeter2Meter();
			cout << it->first << "	" << it->second.x << "	" << it->second.y << "	" << it->second.type_of_point << endl;
			it->second.Meter2Centimeter();

		}
		cout << endl;

	}
	
	HEADING("Covariance of Unknowns");

	IOFormat CleanFmt(2, 0, " ", "\n", "", "");
	cout << "Cx (cm^2) " << "(" << u << " x " << u << "):\n\n";
	cout << Cx.format(CleanFmt) << endl << endl;


	HEADING("Error Ellipses of Unknowns at 39.4%");
	cout << "Point	Semi-Major (cm)	Semi-Minor (cm)	Phi-max-from-N (deg)	Phi-min-from-N (deg)\n\n";

	for (it = pts.begin(); it != pts.end(); it++) {

		if (it->second.type_of_point == 'A') {

			cout << fixed << setprecision(4);
			cout << it->first << "	" << it->second.semi_major << "		" << it->second.semi_minor << "		" << it->second.phi_max << "			" << it->second.phi_min << endl;

		}
	}
	
	cout << endl;

	HEADING("Final Residuals");

	IOFormat CleanFmt2(3, 0, " ", "\n", "", " cm");
	cout << "r " << "(" << n << " x 1" << "):\n\n";
	cout << r.block(0, 0, num_hd_obs, 1).format(CleanFmt2) << endl;

	for (int i = num_hd_obs; i < (num_ha_obs + num_hd_obs); i++) {

		cout << std::right << setw(6) << setprecision(2) << r(i, 0);
		cout << " \"" << endl;
	}
		
	HEADING("Covariance of Residuals");

	IOFormat CleanFmt3(2, 0, " ", "\n", "", "");
	cout << "Cv (row-col 1-6 [units cm^2]) (row-col 7-18 [units (\")^2])	" << "(" << n << " x " << n << "):\n\n";
	cout << Cv.format(CleanFmt3) << endl << endl;
	
	HEADING("Summary of Statistical Summary");

	cout << "		" << "Number of Iterations" << "	=	" << iterations << endl << endl;

	cout << "		" << "Number of Stations" << "	=	" << num_control_points + num_unknown_points << endl << endl;

	cout << "		" << "Number of Observations		" << "	=	" << num_ha_obs + num_hd_obs + num_va_obs << endl;
	cout << "		" << "Number of Unknowns		" << "	=	" << num_unknown_points * type_of_network << endl;
	cout << "		" << "Number of Redundant Observations" << "	=	" << num_ha_obs + num_hd_obs + num_va_obs - (num_unknown_points * type_of_network) << endl << endl;

	cout << "		" << "The a priori variance factor (unitless)	" << "	=	" << apriori << endl;
	cout << "		" << "The a posteriori variance factor (unitless)" << "	=	" << aposteriori << endl << endl;


}