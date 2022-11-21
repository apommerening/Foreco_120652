/*
 *  findNeighbours.cpp
 *  
 *
 *  Created by Arne Pommerening on 20/01/2013.
 *  Copyright 2013 Philodendron International. All rights reserved.
 *
 */

#include <Rcpp.h>
#include <cmath>
#include <vector>

using namespace Rcpp;
using namespace std;

/**
 * Returns translated x2. Modified last on 28.02.2007.
 */
double translateX(double xmax, double x1, double y1, double x2, double y2) {
	double xx2, dx, dy, r, w;
	dx = x2 - x1;
	dy = y2 - y1;
	r = dx * dx + dy * dy;
	xx2 = x2;
	if(dx > xmax * 0.5) {
		w = pow(dx - xmax, 2);
		if(w + pow(dy, 2) < r)
			xx2 = x2 - xmax;
	}
	if(dx < -xmax * 0.5) {
		w = pow(dx + xmax, 2);
		if(w + pow(dy, 2) < r)
			xx2 = x2 + xmax;
	}
	return xx2;
}

/**
 * Returns translated y2. Modified last on 28.02.2007.
 */
double translateY(double ymax, double x1, double y1, double x2, double y2) {
	double yy2, dx, dy, r, w;
	dx = x2 - x1;
	dy = y2 - y1;
	r = dx * dx + dy * dy;
	yy2 = y2;
	if(dy > ymax * 0.5) {
		w = pow(dy - ymax, 2);
		if(pow(dx, 2) + w < r)
			yy2 = y2 - ymax;
	}
	if(dy < -ymax * 0.5) {
		w = pow(dy + ymax, 2);
		if(pow(dx, 2) + w < r)
			yy2 = y2 + ymax;
	}
	return yy2;
}

/**
 * Calculates Euclidean distance between two points with/without implicit translation
 * edge correction. Last modified on 19.02.2007.
 */
double  getEuclideanDistance(int edge, double x1, double y1, double x2, double y2, double xmax, double ymax) {
	double dx, dy, dz, xx2, yy2;
	xx2 = x2;
	yy2 = y2;
	if(edge == 1) { // translation
		xx2 = translateX(xmax, x1, y1, x2, y2);
		yy2 = translateY(ymax, x1, y1, x2, y2);
	}
	dx = xx2 - x1;
	dy = yy2 - y1;
	dz = dx * dx + dy * dy;
	return sqrt(dz);
} 

/**
 * Algorithm to compute the nearest neighbour relations.
 */
// [[Rcpp::export]]
DataFrame findNeighbours(double xmax, double ymax, NumericVector x, NumericVector y, int mi, int edge) {
	double d = 0;
	bool abort;
	int dummy = 0;
	int k = 0;
	/* Number of required neighbours. */
	dummy = mi;
    int na = x.size();
	NumericMatrix distance(na, mi), neighbour(na, mi);
	for(int i = 0; i < na; i++)
		for(int j = 0; j < dummy; j++) {
			distance(i, j) = 1.7E20;
        	neighbour(i, j) = 0;
		}
	for(int i = 0; i < na - 1; i++) {
		for(int j = i + 1; j < na; j++) {
			d = getEuclideanDistance(edge, x[i], y[i], x[j], y[j], xmax, ymax);
			abort = false;
			k = dummy - 1;
			while(abort == false)
				if (k == -1)
					abort = true;
				else if (d < distance(i, k)) 
					k--;
				else abort = true;
			if(k < dummy - 1) {
				for(int l = dummy - 1;l > k + 1; l--) {
					distance(i, l) = distance(i, l - 1);
					neighbour(i, l) = neighbour(i, l - 1);
			    }
				distance(i, k + 1) = d;
				neighbour(i, k + 1) = j;
			}
			abort = false;
			k = dummy-1;
			while (abort == false)
				if(k == -1)
					abort = true;
				else if(d < distance(j, k)) 
					k--;
				else abort = true;
			if (k < dummy - 1) {
				for (int l = dummy - 1; l > k + 1; l--) {
					distance(j, l) = distance(j, l - 1);
					neighbour(j, l) = neighbour(j, l - 1);
				}
				distance(j, k + 1) = d;
				neighbour(j, k + 1) = i;
			}
		}
	}
	return DataFrame::create( Named("distance") = distance, Named ("neighbour") = neighbour);
}


/**
 * Calculates nearest distance of tree to the plot boundary.
 */
double distanceEdge(double xmax, double ymax, double x, double y) {
	double distxy_edge1 = min(xmax - x, ymax - y);
	double distxy_edge2 = min(x - 0, y - 0);
	double distxy_edge = min(distxy_edge1, distxy_edge2);
	return distxy_edge;
}

/**
 * Distance to neighbour must be smaller than distance to the edge.
 */
int neighbourEdgeOk(double xmax, double ymax, double x, double y, double distance) {
	double di = distanceEdge(xmax, ymax, x, y);
	int ateb = 0;
	if(di > distance)
		ateb = 1;
    return ateb;
}

/**
 * Calculating Horvitz-Thompson weights.
 */
double calcRepFactor(int edgeCorr, double bufferWidth, double xmax, double ymax, double x, double y, double distance) {
	double ci = 0;
	double RF = 1;
	if(edgeCorr == 3) {
		double di = distanceEdge(xmax, ymax, x, y);
		if(di <= bufferWidth)
			RF = 0;
	}
	if(edgeCorr > 3) {
	  if(edgeCorr == 4)
	    ci = distance;
	  double area = ((xmax - 2 * ci) * (ymax - 2 * ci)) / 10000;
	  RF = neighbourEdgeOk(xmax, ymax, x, y, distance) * 1 / area;
	}
	return RF;
}

// [[Rcpp::export]]
DataFrame calcRepFactors(int edgeCorr, double bufferWidth, double xmax, double ymax, NumericVector x, NumericVector y, NumericVector dist1, NumericVector dist4) {
	int n = x.size();
    NumericVector rf1(n), rf4(n);
	for(int i = 0; i < n; i++) {
	  rf1[i] = calcRepFactor(edgeCorr, bufferWidth, xmax, ymax, x[i], y[i], dist1[i]);
  	  rf4[i] = calcRepFactor(edgeCorr, bufferWidth, xmax, ymax, x[i], y[i], dist4[i]);
	}
	return DataFrame::create( Named("rf1") = rf1, Named ("rf4") = rf4);
}
