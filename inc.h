/*
 * inc.h
 *
 *  Created on: 21/05/2014
 *      Author: yan
 */

#ifndef INC_H_
#define INC_H_

#include <iostream>
//#include <vector>
//#include <math.h>
//#include <complex>
//#include <fftw3.h>
//#include <unistd.h>
//#include <stdlib.h>
//#include <string.h>
#include <fstream>
#include <valarray>
//#include <gsl/gsl_multifit.h>
//#include <stdbool.h>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/contrib/contrib.hpp>

#define PLOT 0

using namespace std;
using namespace cv;
//#define GNUPLOT "gnuplot -persist"
static const double pi = 3.141592653589793238460;
static const double C = 3e+8; // velocity of light in m/sec
double PRF = 10.0e3; //pulse repetition frequency
const complex<double> J(0.0, 1.0);


#endif /* INC_H_ */
