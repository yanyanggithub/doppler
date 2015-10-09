/*
 * util.h
 *
 *  Created on: 22/05/2014
 *      Author: yan
 */

#ifndef UTIL_H_
#define UTIL_H_
#include "inc.h"



//Normalize to [-180,180):
inline double constrainAngle(double x){
	x = fmod(x + pi,2*pi);
	if (x < 0)
		x += 2*pi;
	return x - pi;
}

// convert to [-360,360]
inline double angleConv(double angle){
	return fmod(constrainAngle(angle),2*pi);
}

inline double angleDiff(double a,double b){
	double dif = fmod(b - a + pi,2*pi);
	if (dif < 0)
		dif += 2*pi;
	return dif - pi;
}

inline double unwrap(double previousAngle,double newAngle){
	return previousAngle - angleDiff(newAngle,angleConv(previousAngle));
}

void unwrap(vector<double> &x)
{
	if(x.size() > 1)
	{
		for(size_t i = 1; i < x.size(); i++)
		{
			x[i] = unwrap(x[i-1], x[i]);
		}
	}
}

vector<double> sampling(double min, double interval, double max)
{
	vector<double> samples;
	int N = floor((max - min)/interval);
	double t = 0.0;
	for(int i = 0; i <= N; i++)
	{
		t = min + interval * i;
		samples.push_back(t);
	}

	return samples;
}

void quadraticfit(vector<double> &dx, vector<double> &dy, vector<double> &coeffs)
/* least square estimation to fit a quadratic function */
{
	double S00, S01, S10, S11, S20, S21, S30, S40 = 0.0;

	S00 = (dx.size());
	for(size_t i = 0; i < dx.size(); i++)
	{
		double dxi = dx[i];
		double dyi = dy[i];
		S10 += dxi;
		S20 += pow(dxi, 2);
		S30 += pow(dxi, 3);
		S40 += pow(dxi, 4);

		S01 += dyi;
		S11 += dxi*dyi;
		S21 += pow(dxi, 2)*dyi;
	}

	double a, b, c = 0.0;
	a = (S01*S10*S30 - S11*S00*S30 - S01* pow(S20,2)
		       + S11*S10*S20 + S21*S00*S20 - S21*pow(S10,2))
		    /(S00*S20*S40 - pow(S10,2)*S40 - S00*pow(S30,2) + 2*S10*S20*S30 - pow(S20,3));
	b = (S11*S00*S40 - S01*S10*S40 + S01*S20*S30
		       - S21*S00*S30 - S11*pow(S20,2) + S21*S10*S20)
		    /(S00*S20*S40 - pow(S10,2)*S40 - S00*pow(S30,2) + 2*S10*S20*S30 - pow(S20,3));
	c = (S01*S20*S40 - S11*S10*S40 - S01*pow(S30,2)
		       + S11*S20*S30 + S21*S10*S30 - S21*pow(S20,2))
		    /(S00*S20*S40 - pow(S10,2)*S40 - S00*pow(S30,2) + 2*S10*S20*S30 - pow(S20,3));

	coeffs.push_back(a);
	coeffs.push_back(b);
	coeffs.push_back(c);
}

double polyval(vector<double> q, double x)
{
	long double val = 0.0;
	int degree = q.size();
	for(int i = 0; i < degree; i++)
	{
		val += pow(x, degree-i-1) * q[i];
	}

	return val;
}

double gaussian (double x, double mu, double sigma) {
  return exp( -(((x-mu)/(sigma))*((x-mu)/(sigma)))/2.0 );
}

int nextpow2(double x)
{
	int p = 0;
	while(pow(2.0,p) <= x)
	{
		p = p + 1;
	}

	return p;
}

void initVector(int rows, int cols, vector<valarray<complex<double> > > &y)
{
	y.clear();
	complex<double> zero(0.0,0.0);
	for(int i = 0; i < rows; i++)
	{
		valarray<complex<double> > yi(cols);
 		y.push_back(yi);
	}
}

void initVector(int N, vector<complex<double> > &y)
{
	y.clear();
	complex<double> zero(0.0,0.0);
	for(int i = 0; i < N; i++)
	{
		y.push_back(zero);
	}
}

bool IsPowerofTwo(int x)
{
	return (x != 0) && ((x & (x - 1)) == 0);
}

template <typename T>
void transpose(vector< valarray<T> > &x, vector< valarray<T> > &xt )
{

	int M = x.size();
	int N = x[0].size();
	initVector(N, M , xt);
	for(int i = 0; i < M; i++)
	{
		for(int j = 0; j < N; j++)
		{
			xt[j][i] = x[i][j];
		}
	}
}

template <typename T>
void fft(valarray<T> &x, bool radix)
/* fast Fourier transform for a vector
 * time domain -> frequency domain */
{
	const size_t N = x.size();

	if(radix)
    /* if signal size is the power of 2, use cooley-turkey*/
	{
		if (N <= 1) return;
		// divide
		valarray<T> even = x[std::slice(0, N/2, 2)];
		valarray<T>  odd = x[std::slice(1, N/2, 2)];

		// conquer
		fft(even, radix);
		fft(odd, radix);

		// combine
		for (size_t k = 0; k < N/2; ++k)
		{
			T t = std::polar(1.0,-2.0*pi*(double)k/(double)N)*odd[k];
			x[k    ] = even[k] + t;
			x[k+N/2] = even[k] - t;
		}
	} else
	/* slow ft */
	{
		valarray<T> out(N);
		for(size_t k = 0; k < N; k++)
		{
			T Xk(0.0, 0.0);
			for(size_t j = 0; j < N; j++)
			{
				double img = -2.0 * pi * (double)j*(double)k/(double)N;
				T v(0, img);
				Xk = Xk + x[j] * exp(v);
			}
			out[k] = Xk;
		}
		x = out;
	}

}

template <typename T>
void ifft(valarray<T> &x, bool radix)
/* inverse Fourier transform
 * frequency domain -> time domain */
{
	if(radix)
	/* if signal size is the power of 2, use cooley-turkey*/
	{
		// conjugate the complex numbers
		x = x.apply(std::conj);

		// forward fft
		fft(x, radix);

		// conjugate the complex numbers again
		x = x.apply(std::conj);

		// scale the numbers
		x /= (double)x.size();
	}else
	{
		/* slow ift */
		int N = x.size();
		valarray<T> out(N);
		for(int j = 0; j< N; j++)
		{
			T xj(0.0, 0.0);
			for(int k = 0; k < N; k++)
			{
				double img = 2.0 * pi * (double)j*(double)k/(double)N;
				T v(0, img);
				xj = xj + (1.0/(double)N) * x[k] * exp(v);
			}
			out[j] = xj;
		}
		x = out;
	}
}

template <typename T>
void conv(valarray<T> &h, valarray<T> &y, valarray<T> &out)
{
	int M = y.size();
	int L = h.size();
	for(int i = 0; i < M+L; i++)
	{
		int i1 = i;
		for(int k = 0; k < L; k++)
		{
			if(i1 >= 0 && i1 < M)
				out[i] = out[i] + (y[i1]*h[k]);
			i1 = i1-1;
		}
	}
}

template <typename T>
void conv(valarray<T> &h,
		vector< valarray<T> > &y,
		vector< valarray<T> > &out, bool rows)
{
	int M = y.size();
	int N = y[0].size();
	int L = h.size();
	if(rows)
	/*convolution on rows*/
	{
		for(int i = 0; i < M; i++)
		{
			conv(h, y[i], out[i]);
		}
	}
	else
	/*convolution on columns*/
	{
		for(int j = 0; j < N; j++)
		{
			for(int i = 0; i < M+L; i++)
			{
				int i1 = i;
				for(int k = 0; k < L; k++)
				{
					if(i1 >= 0 && i1 < M)
						out[i][j] = out[i][j] + (y[i1][j]*h[k]);
					i1 = i1-1;
				}
			}
		}
	}
}


double rand_uniform()
{
	return (double) rand()/(double)RAND_MAX;
}

template <typename T>
void git_rotate(vector<valarray<T> > &x, int num_places)
{
	int M = x.size();
	int N = x[0].size();

	if(N > 1)
    /* rotate columns */
	{
		num_places = num_places%N;
		for(int i = 0 ; i < M; i++)
		{
			for(int j = 0 ; j < N-num_places; j++)
			{
				T temp = x[i][j];
				x[i][j] = x[i][N-num_places+j];
				x[i][N-num_places+j] = temp;
			}
		}
	}
	/* rotate rows */
	else if(M > 1)
	{
		num_places = num_places%M;
		for(int i = 0 ; i < M-num_places; i++)
		{
			for(int j = 0 ; j < N; j++)
			{
				T temp = x[i][j];
				x[i][j] = x[M-num_places+i][j];
				x[M-num_places+i][j] = temp;
			}
		}
	}
}

void peakinterp(double &amp, double &del_index,
		double &z1, double &z2, double &z3)
/* peakinterp performs a quadratic interpolation of the peak
 * defined by three consecutive data values in the three-element vector z.
 * the middle value z2 must be the largest element.
 * peakinterp returns the interpolated peak amplitude and the interpolated
 * peak location relative to the center element in the range bins*/
{
	del_index = -0.5 * ((z3 - z1))/(z1 - 2*z2 + z3);
	int k = del_index;
	amp = 0.5 * ( (k-1)*k*z1 - 2*(k-1)*(k+1)*z2 + (k+1)*k*z3 );
}

double db(complex<double> x, string label)
/* converts energy or power measurement to decibels */
{
	double val = 0.0;
	double ab = abs(x);
	if(label == "voltage")
	{
		val = (ab == 0) ? 0 : 20*log10(ab);
	} else if(label == "power")
	{
		val = (ab == 0) ? 0 : 10*log10(ab);
	} else {
		cout << "invalid label of db" << endl;
	}

	return val;
}

vector<double> hamming(int L)
/* L-point symmetric Hamming window in the column vector w*/
{
	vector<double> w;
	for(int i = 0; i < L; i++)
	{
		double val = 0.54 - 0.46 * cos(2.0*pi*(double)i/(double)(L-1));
		w.push_back(val);
	}
	return w;
}


#endif /* UTIL_H_ */
