/*
 * radar.h
 *
 *  Created on: 21/05/2014
 *      Author: yan
 */

#ifndef RADAR_H_
#define RADAR_H_
#include "inc.h"
#include "util.h"

valarray<complex<double> > chirp_signal(double T, double W, double p)
/*
 * generate a sampled chirp signal
 */
{

	double delta_t = 1.0 / (p*W);
	int N = p*T*W;
	valarray<complex<double> > x(N);
	for(int i = 0; i < N; i++)
	{
		double d = pow((delta_t*i-(double)(N-1)/2.0/p/W), 2);
		double j = pi* W/T *d;
		complex<double> t(0,j);
//		x.push_back(expt);
		x[i] = exp(t);
	}

	return x;
}

class radar
/* Make pulse Doppler radar data
 * + velocity in meters/sec, distance in m, times in sec, BW in Hz
 * + assumes each pulse is constant (complex) amplitude
 * + will accommodate up to quadratic phase pulses
 * + allows distributed targets
 */
{
public:
	valarray<complex<double> > x; //baseband single pulse waveform
	double fs; // sampling frequency of input pulse, Hz
	vector<double> T_0; // start time(s) of input pulse(s), sec
						// number of pulses in burst assumed = length(g)
	vector<double> g; //complex gain(s) of pulse(s)
	vector<double> T_out; // 2-vector [T_min, T_max] defines output window delay
						 // times w.r.t start of pulse
	double T_ref; // system "reference" time
	double fc; // center frequency of the radar
	vector<double> ranges; // vector of ranges to target(s), meters
	vector<double> snr; // vector of target SNRs (unit noise power assumed)
	vector<double> v; // vector of target velocities, m/sec

	void setup(valarray<complex<double> > x, double fs, vector<double> T_0,
			vector<double> g, vector<double> T_out,	double T_ref,
			double fc,	vector<double> ranges, vector<double> snr,
			vector<double> v);
	vector< valarray<complex<double> > > simulate();
	void add_noise(vector< valarray<complex<double> > > &y);
};

void radar::setup(valarray<complex<double> > x, double fs, vector<double> T_0,
		vector<double> g, vector<double> T_out,	double T_ref,
		double fc,	vector<double> ranges, vector<double> snr,
		vector<double> v)
{
	this->x = x;
	this->fs = fs;
	this->T_0 = T_0;
	this->g = g;
	this->T_out = T_out;
	this->T_ref = T_ref;
	this->fc = fc;
	this->ranges = ranges;
	this->snr = snr;
	this->v = v;
}

vector< valarray<complex<double> > > radar::simulate()
{
	vector< valarray<complex<double> > > y;
	int Mx = x.size();
	double delta_t = 1.0/fs; // sampling interval (sec)

	//output sampling times (sec)
	vector<double> t_y = sampling(T_out[0], delta_t, T_out[1]);
	double T_p = Mx * delta_t; // length of input pulse (sec)

	/* determine the quadratic phase modulation parameters for later
	 * interpolation of pulse samples*/
	vector<double> t_x;
	for(int i = 0; i < Mx; i++)
	{
		t_x.push_back(i*delta_t*1e9); // use msec for quadratic fitting
	}
	vector<double> x_ph;
	for(int i = 0; i < Mx; i++)
	{
		x_ph.push_back(arg(x[i]));
	}
	unwrap(x_ph);
	vector<double> q; // coefficients of the quadratic function
	quadraticfit(t_x, x_ph, q);

	int Mr = t_y.size(); //rows of output matrix y
	int Nr = g.size(); // cols of output matrix y
	initVector(Mr, Nr, y);

	/*'i' loops over the number of targets*/
	for(size_t i = 0; i < ranges.size();i++)
	{
		double ri = ranges[i];
		double vi = v[i];
		double f_doppler = 2.0 * vi * fc/C;


		/*'j' loops over the number of pulses*/
		for(size_t j = 0; j < g.size(); j++)
		{
			double r_at_T_0 = ri - vi*T_0[j];

			/**
			 * compute start and end time of reflected pulse at receiver
			 * ensure that it falls at least partially within the range window
			 */
			double tau = 2.0 * r_at_T_0 / (C + vi);
			double tmax = tau + T_p;
//			cout << tau << " " << tmax << endl;

			if(tau >= T_out[1] || tmax <= T_out[0])
			{
				cout << "Echo from target " << i << "at range " << ri << "km" << endl;
				cout << "is COMPLETELY OUT OF the range window" << endl;
				cout << "on pulse " << j << endl;
			} else{
				/* place scaled, range-delayed, Doppler shifted pulse
				 * into output matrix
				 * unit noise power and unit nominal pulse amplitude assumed
				 * to get amplitude from SNR*/
				double amp = pow(10, snr[i]/20.0);
				for(size_t t = 0; t < t_y.size(); t++)
				{
					/*
					 * figure out which sample locations in the output grid contain
					 * reflected pulse
					 * */
					double t_vals_t = t_y[t] - tau;
					if(t_vals_t >= 0 && t_vals_t < T_p)
					{
						complex<double> t1(0, -2.0*pi*fc*tau);
						complex<double> t2(0, 2.0*pi*f_doppler*t_y[t]);
						// convert t_vals_t to msec for polyval
						complex<double> t3(0, polyval(q, t_vals_t*1e9));
//						cout << exp(t1) << " " << exp(t2) << " " << exp(t3) << endl;
						y[t][j] = y[t][j] + (amp * g[j] * exp(t1)) * exp(t2) * exp(t3);
					}
				}
			}
		}
	}

	add_noise(y);
	return y;
}

void radar::add_noise(vector< valarray<complex<double> > > &y)
{
	int My = y.size();
	int Ny = y[0].size();

//	srand( time( NULL ) );
	srand(0);
	Mat real(My, Ny, CV_64FC1);
	Mat imag(My, Ny, CV_64FC1);
	Mat r(My, Ny, CV_64FC1);
	randn(real, Scalar(0), Scalar(1));
	randn(imag, Scalar(0), Scalar(1));
	randn(r, Scalar(0), Scalar(1));

	/* adding noise */
	for(int i = 0; i< My; i++)
	{
		for(int j = 0; j < Ny; j++)
		{
			complex<double> white_noise(real.at<double>(i, j), imag.at<double>(i, j));
			y[i][j] += y[i][j] + (1.0/(double)sqrt(2.0)) * white_noise;
		}
	}

	/* creating clutter
	 * create log-normal (ground) "clutter" with specified C/N
	 * and log-normal standard deviation for amplitude, uniform phase
	 * clutter is uncorrelated in range, fully correlated in pulse
	 * */
	int CN = 20; // clutter-to-noise ratio in first bin (dB)
	int SDxdB = 3; // in dB (is not the sigma of the complete clutter)
	vector< valarray<complex<double> > > ncc;
	for(int i = 0; i < My; i++)
	{
		valarray<complex<double> > ncci(Ny);
		for(int j = 0; j < Ny; j++)
		{
			double val = pow(10.0, SDxdB * r.at<double>(i,j)/10.0);
			complex<double> v1(0, 2.0*pi*rand_uniform());
			ncci[j] = val * exp(v1);
		}
		ncc.push_back(ncci);
	}

	/* correlating and adding clutter
	 * force the power spectrum shape to be Gaussian*/
	vector<double> G;
	double gau[5];
	for(int i = 0; i < 5; i++)
	{
		gau[i] =  exp(- (i*i)/1.0);
	}

	for(int i = 0; i < Ny; i++)
	{
		if(i < 5)
			G.push_back(gau[i]);
		else if(i >= Ny-5+1)
		{
			G.push_back(gau[Ny-i]);
		} else
			G.push_back(0.0);
	}

	complex<double> mean(0.0,0.0);
	complex<double> var(0.0,0.0);
	for(int i = 0; i < My; i++)
	{
		fft(ncc[i], IsPowerofTwo(Ny));

		for(int j = 0 ; j < Ny; j++)
		{
			ncc[i][j] = ncc[i][j] * G[j];
		}

		ifft(ncc[i], IsPowerofTwo(Ny));
		mean += ncc[i].sum();
	}
	mean = mean/(double)(My*Ny);
//	cout << mean << endl;
	for(int i = 0; i < My ; i++)
	{
		for(int j = 0; j < Ny; j++)
		{
			complex<double> nc = ncc[i][j];
			var += (nc-mean) * conj(nc-mean);
		}
	}

	var = var/(double)(My*Ny);
//	cout << var << endl;

	/* rescale clutter to have desired C/N ratio */
	double scale = sqrt( (double)pow(10, CN/10)/var.real());
//	cout << scale << endl;
	for(int i = 0; i < My ; i++)
	{
		/* weight the clutter power in range for assume R^2 (beam-limited) loss */
		double cweight = T_out[0] * pow( T_out[0] + i * (1.0/fs) , -1);
		for(int j = 0; j < Ny; j++)
		{
			complex<double> nc = ncc[i][j];
			y[i][j] = y[i][j] + cweight * scale * nc;
		}
	}
}

#endif /* RADAR_H_ */
