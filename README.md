# Geodetic network adjustment (C++)

Parametric least-squares adjustment of a 2D survey network from horizontal angle and
distance observations. Written in C++ with Eigen.

Given fixed control points, approximate coordinates for the unknown points and the
observations with their standard deviations, the program linearizes the observation
equations, iterates the adjustment to convergence, and reports:

- adjusted coordinates of the unknown points
- variance-covariance matrices of the unknowns, observations and residuals
- residuals for every observation
- absolute error ellipses for each unknown point
- the a posteriori variance factor

## Results on the sample network

Three control points (A, B, C), two unknown points (1, 2), 12 horizontal angles and
6 distances (18 observations, 4 unknowns, 14 redundant). Converged in 4 iterations.

| Point | X (m) | Y (m) | Semi-major (cm) | Semi-minor (cm) |
|---|---|---|---|---|
| 1 | 5656.9414 | 2475.5609 | 3.04 | 2.74 |
| 2 | 663.8917 | 2944.0507 | 5.42 | 4.84 |

Error ellipses are at 39.4% confidence (one standard deviation in 2D).

The a posteriori variance factor came out at 29.44 against an a priori value of 2.25,
a ratio of about 13. That points to the stated observation precisions (1.5" for
angles, 2 cm for distances) being too optimistic for this data set.

## Method

1. Build the observation covariance matrix from the stated precisions.
2. Form the design matrix (partial derivatives of each angle and distance equation with
   respect to the unknown coordinates) and the misclosure vector at the current estimate.
3. Solve the normal equations, update the coordinates and repeat until the corrections
   fall below a threshold.
4. Propagate the covariance of the unknowns and compute error ellipses from its
   eigenvalues.

The network classes derive from an abstract `PositioningNetwork` interface, so other
network types (levelling, 3D, image observations) can be added with the same structure.

## Build and run

Requires a C++17 compiler, CMake 3.10+ and [Eigen 3](https://eigen.tuxfamily.org).

```bash
cmake -S . -B build              # add -DEIGEN3_INCLUDE_DIR=/path/to/eigen if CMake can't find it
cmake --build build
./build/gna data/apriori.txt data/observations.txt
```

Run it from the repository root so the default data paths resolve.

## Input format

`data/apriori.txt`: the first line `! 2` marks a 2D network. Then one point per line:

```
C  <id>  <x>  <y>     control (fixed) point, metres
A  <id>  <x>  <y>     approximate coordinates of an unknown point, metres
```

`data/observations.txt`: one observation per line:

```
H  <from>  <at>  <to>  <deg> <min> <sec>  <std, arcsec>    horizontal angle
D  <from>  <to>  <distance, m>  <std, m>                   horizontal distance
```

## Layout

```
src/        main.cpp (driver), GeoNetwork.cpp (adjustment)
include/    GeoNetwork.h, PositioningNetwork.h (abstract interface)
data/       sample network
```

Written in 2019 for a data analysis course in Geomatics Engineering at the University of Calgary.
