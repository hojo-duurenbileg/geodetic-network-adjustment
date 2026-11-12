#ifndef POSITIONINGNETWORK_H
#define POSITIONINGNETWORK_H

//**
// Summary:
//    An abstract network class to layout the general concepts from which more
//	  specific network classes can be derived. Provides the basic functionality/interface.
//
// History:
//    Feb 10/18 - Created by Khongorzul Duurenbileg
//

#include <string>

class PositioningNetwork
{
public:
	// mutators

	//**
	// Summary:
	//    Expression for setting various types of points (1D, 2D, 3D, interest, fiducial, control...).
	//
	// Arguments:
	//    
	// Returns:
	//    
	virtual void SetPoints(std::string filename) = 0;

	//**
	// Summary:
	//    Expression for setting various types of observables (angles, distances, bearings, image observations...).
	//
	// Arguments:
	//    
	// Returns:
	//    
	virtual void SetObservables(std::string filename) = 0;
	//**
	// Summary:
	//    Expression for setting various types of observables (angles, distances, bearings, image observations...).
	//
	// Arguments:
	//    
	// Returns:
	//    
	virtual void ReadData(std::string filename) = 0;

private:



};


#endif // !POSITIONINGNETWORK_H
