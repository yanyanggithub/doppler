/*
 * simulator.cpp
 *
 *  Created on: 21/05/2014
 *      Author: yan
 */

#ifndef SIMULATOR_CPP_
#define SIMULATOR_CPP_
#include "inc.h"
#include "util.h"
#include "radar.h"

void write_signal(valarray<complex<double> > &s, vector<double> &sx, string filename)
{
	ofstream file(filename.c_str());
	file << "#t real imag" << endl;
	for(size_t i=0; i< s.size(); i++){
		file << sx[i] << ' ' << s[i].real() << ' ' << s[i].imag() << endl;
	}
	file.close();
}

void write_data_plot(vector< vector<double> > &data, vector<double> &x_axis,
		string filename)
{
	ofstream file(filename.c_str());

	int Mr = data.size();
	int Nr = data[0].size();
	for(int i = 0; i < Mr; i++ )
	{
		file << x_axis[i] << ' ';
		for(int j = 0; j < Nr; j++)
		{
			file << data[i][j] << ' ';
		}
		file << endl;
	}
	file.close();
}

void write_3d_data_plot(vector< vector<double> > &data, vector<double> &x_axis,
		string filename, bool columndata)
{
	ofstream file(filename.c_str());

	int Mr = data.size();
	int Nr = data[0].size();
	if(columndata)
	{
		for(int i = 0; i < Nr; i++ )
		{
			for(int j = 0; j < Mr; j++)
			{
				file << i << ' ' << x_axis[j] << ' ' << data[j][i] << '\n';
			}
			file << endl;
		}
	} else
	{
		for(int i = 0; i < Mr; i++ )
		{
			for(int j = 0; j < Nr; j++)
			{
				file << i << ' ' << x_axis[j] << ' ' << data[i][j] << '\n';
			}
			file << endl;
		}
	}
	file.close();
}

void get_a_radar(radar &a_radar)
{
	/* form radar chirp pulse */
	double T = 10e-6; // pulse length, sec
	double W = 10e6; // chirp bandwidth, Hz
	double fs = 12e6; // chirp sampling rate, Hz
	valarray<complex<double> > s = chirp_signal(T, W, fs/W);

	int Np = 20; // 20 pulses
	double PRI = 1.0/PRF; // PRI in sec
	vector<double> T_0; // relative start times of pulse, in sec
	vector<double> g; // gains of pulses
	for(int i = 0; i < Np; i++)
	{
		T_0.push_back(PRI * (double)i);
		g.push_back(1.0);
	}
	vector<double> T_out; // start and end times of range window in sec
	T_out.push_back(12e-6);
	T_out.push_back(40e-6);
	double T_ref = 0.0; //system reference time in usec
	double fc = 10e9; // RF frequency in Hz; 10 GHz is x-band

	/* compute unambiguous Doppler interval in m/sec
	 * compute unambiguous range interval in meters
	 */
	double vua = C*PRF/(2*fc);
	double rmin = C*T_out[0]/2;
	double rmax = C*T_out[1]/2;
	double rua = C/2/PRF;
//	int Ntargets = 4;
//	double del_R = (C/2.0) * (1.0/fs)/1e3;
	cout << "the unambiguous velocity interval is " << vua << "m/s" << endl;
	cout << "the range window starts at " << rmin/1e3 << "km" << endl;
	cout << "the range window ends at " << rmax/1e3 << "km" << endl;
	cout << "the unambiguous range interval is " << rua/1e3 << "km" << endl;

	vector<double> ranges, SNR, vels;
	ranges.push_back(2e3); SNR.push_back(-3); vels.push_back(-0.4*vua); //target 1
	ranges.push_back(3.8e3); SNR.push_back(5); vels.push_back(-0.2*vua); //target 2
	ranges.push_back(4.4e3); SNR.push_back(10); vels.push_back(0.2*vua); //target 3
	ranges.push_back(4.4e3); SNR.push_back(7); vels.push_back(0.4*vua); //target 4

	/* compute relative RCS using the idea that SNR is proportional to
	 * RCS/R^4
	 */
	vector<double> rel_RCS;
	double max_rcs = 0.0;
	for(size_t i = 0; i < SNR.size(); i++)
	{
		double val = pow(10, SNR[i]/10) * pow(ranges[i], 4);
		if(val > max_rcs)
		{
			max_rcs = val;
		}
		rel_RCS.push_back(val);
	}
	// convert to power
	for(size_t i = 0; i < rel_RCS.size(); i++)
	{
		double temp = rel_RCS[i];
		rel_RCS[i] = db(temp/max_rcs, "power");
	}

	a_radar.setup(s, fs, T_0, g, T_out, T_ref, fc, ranges, SNR, vels);
}

void get_radar_signal(radar &a_radar, vector< valarray<complex<double> > > &y)
{
	y = a_radar.simulate();

	/* output the data for plots */
	vector<double> x_axis;
	for(size_t i = 0; i < y.size(); i++)
	{
		double val = ((double)C/2.0) * ((i * (1.0/a_radar.fs)) + a_radar.T_out[0]) / (double)1e3;
		x_axis.push_back(val);
	}

	int My = y.size();
	int Ny = y[0].size();
	vector< vector<double> > ydB, ydBscaled;
	double ydBmax = 0.0;
	for(int i = 0; i < My; i++)
	{
		vector<double> ydbi;
		vector<double> ydbscaledi;
		ydbi.clear();
		ydbscaledi.clear();
		for(int j = 0; j < Ny; j++)
		{
			ydbi.push_back(db(y[i][j], "voltage"));
			double val = abs(y[i][j]);
			if(val > ydBmax)
				ydBmax = val;
			ydbscaledi.push_back(val);
		}
		ydB.push_back(ydbi);
		ydBscaled.push_back(ydbscaledi);
	}
	write_data_plot(ydB, x_axis, "y.dat");

	for(int i = 0; i < My; i++)
	{
		for(int j = 0; j < Ny; j++)
		{
			ydBscaled[i][j] = db(ydBscaled[i][j]/ydBmax, "voltage");
		}
	}
	write_3d_data_plot(ydBscaled, x_axis, "ydb.dat", true);
}

void doppler_process(vector<valarray<complex<double> > > &y,
		vector<valarray<complex<double> > > &yp, vector<double> &w,  int Lfft)
{
//	vector<valarray<complex<double> > > yp;
	int My = y.size();
	int Ny = y[0].size();
	int N = (Lfft <= Ny) ? Ny : Lfft;
	for(int i = 0; i < My; i++)
	{
		valarray<complex<double> > in(N);
		vector<double> YdBi;

		YdBi.clear();
		for(int j = 0; j < Ny; j++)
		{
			in[j] = conj(y[i][j]) * w[j];
		}

		fft(in, IsPowerofTwo(N));
		in = in * in.apply(std::conj);
		yp.push_back(in);
	}
	git_rotate(yp, Lfft/2);
}

void range_doppler_plot(vector<valarray<complex<double> > > &y,
		vector<double> &range, int Lfft, string filename )
{
//	int My = y.size();
	int Ny = y[0].size();
	vector<double> w = hamming(Ny);

	vector<valarray<complex<double> > > yp;
	vector< vector<double> > YdB;
	doppler_process(y, yp, w, Lfft);
	int M = yp.size();
	int N = yp[0].size();
	cout << yp.size() << " " << yp[0].size() << endl;

	double maxdb = 0.0;

	for(int i = 0; i < M; i++)
	{
		vector<double> YdBi;

		for(int j = 0; j < N; j++)
		{
			double val = db(yp[i][j], "power");
			if(val > maxdb)
				maxdb = val;
			YdBi.push_back(val);
		}
		YdB.push_back(YdBi);
	}

	cout << filename << ":  " << maxdb << endl;
	write_3d_data_plot(YdB, range, filename, false);
}

void range_plot(vector<valarray<complex<double> > > &y,
		vector<double> &range, int Lfft, string filename )
{
	int M = y.size();
	int N = y[0].size();
	cout << y.size() << " " << y[0].size() << endl;
	vector< vector<double> > YdB;
	double maxdb = 0.0;

	for(int i = 0; i < M; i++)
	{
		vector<double> YdBi;

		for(int j = 0; j < N; j++)
		{
			double val = db(y[i][j], "power");
			if(val > maxdb)
				maxdb = val;
			YdBi.push_back(val);
		}
		YdB.push_back(YdBi);
	}

	cout << filename << ":  " << maxdb << endl;
	write_3d_data_plot(YdB, range, filename, true);
}

void process_radar_signal(radar &a_radar, vector< valarray<complex<double> > > &y)
{
	vector<double> sx;
	for(size_t i = 0; i < a_radar.x.size(); i++)
	{
		double val = ((double)1e6/a_radar.fs)*i;
		sx.push_back(val);
	}
#if PLOT
	write_signal(a_radar.x, sx, "chirp.dat");
#endif

//	double PRI = 1.0/PRF;

	//compute unambiguous Doppler interval in m/sec
	//compute unambiguous Dopper range interval in meters
	double vua = (-1) * C*PRF/(2*a_radar.fc);
//	double rmin = C*a_radar.T_out[0]/2;
//	double rmax = C*a_radar.T_out[1]/2;
//	double rua = C/2/PRF;

	//convert range samples to absolute range units
	int My = y.size();
	int Ny = y[0].size();
	vector<double> range;
	for(int i = 0; i < My; i++)
	{
		double val = (double(C)/2.0)* ((double)i*(1.0/a_radar.fs) + a_radar.T_out[0])/1e3;
		range.push_back(val);
	}
	vector<int> pulses;
	for(int i = 0; i < Ny; i++)
	{
		pulses.push_back(i);
	}

	//force oversize FFT, and compute doppler scale factor
	int p = nextpow2(Ny);
	int Lfft = pow(2.0, p+3);
	vector<double> doppler;
	for(int i = 0; i < Lfft; i++)
	{
		double val = ((double)i/(double)Lfft - 0.5) * vua;
		doppler.push_back(val);
	}

	/* dopper process and square-law detect the whole
	 * unprocessed array
	 * using Hamming window throughout*/
#if PLOT
	range_doppler_plot(y, doppler, Lfft, "doppler_range.dat");
#endif

	/* processing the data
	 * 1. pulse compression
	 * 	  using time-domain hamming weighting of the
	 * 	  impulse response for range sidelobe control */
	int Ls = a_radar.x.size();
	vector<double> hw = hamming(Ls);
	valarray<complex<double> > h1(Ls);
	for(int i = 0; i < Ls; i++)
	{
		h1[i] = conj(a_radar.x[Ls-1-i]) * hw[i];
	}
	int Myp = My + Ls -1;
	int Nyp = Ny;
	vector< valarray<complex<double> > > yp;
	initVector(Myp, Nyp, yp);
	conv(h1, y, yp, false);
	vector<double> rangep;
	for(int i = 0; i < Myp; i++)
	{
		double rval = (C/2.0)*((double(i)-(Ls-1))*(1.0/a_radar.fs) + a_radar.T_out[0])/1e3;
		rangep.push_back(rval);
	}
#if PLOT
	range_doppler_plot(yp, doppler, Lfft, "doppler_range_compressed.dat");
#endif

	/* 2. apply three-pulse canceller in each range bin to raw data */
	valarray<complex<double> > h2(3);
	h2[0] = complex<double>(1, 0); h2[1] = complex<double>(-2, 0);
	h2[2] = complex<double>(1, 0);
	vector< valarray<complex<double> > > ypm;
	int Nyp2 = Nyp + h2.size() - 1;
	initVector(Myp, Nyp2, ypm);
	conv(h2, yp, ypm, true);
	vector< valarray<complex<double> > > YPM;
	vector<double> w = hamming(Ny);
	doppler_process(ypm, YPM, w, Lfft);
#if PLOT
	range_doppler_plot(ypm, doppler, Lfft, "doppler_range_canceller.dat");
#endif

	/*3. search for the range bins with targets
	 *   first noncoherently integrate across the frequency bins*/
	int Mypm = YPM.size();
//	int Nypm = YPM[0].size();
	valarray<double> YPMrange(Mypm);
	vector<double> ranges_sort;
	for(int i = 0; i < Mypm; i++)
	{
		double sumr = YPM[i].sum().real();
		YPMrange[i] = sumr;
		ranges_sort.push_back(sumr);
	}
	sort(ranges_sort.begin(), ranges_sort.end());
	/* median of YPMrange */
	int md_ind = (Mypm%2) ? Mypm/2 : (Mypm/2+1);
	double Nrange = ranges_sort[md_ind];
	/* threshold 8x (9DB) above noise estimate */
	double Trange = 8*Nrange;

	/* loop identifies which range bins have local peaks above
	 * the threshold. It also set up a vector */
	vector<pair<int, double> > spikesr;
	for(int i = 2; i < Mypm-1; i++)
	{
		if( (YPMrange[i] > YPMrange[i+1]) && (YPMrange[i]>YPMrange[i-1])
		             && (YPMrange[i] > Trange) )
		{
			spikesr.push_back(make_pair(i, YPMrange[i]));
		}
	}

	vector<vector<double> > targets;
	/*  4. find the Doppler peak(s) for each range bin having a target(s).  Keep
		adjoining Doppler values as well to support subsequent interpolation
	*/
 	int Mspikesr = spikesr.size();
	for(int i = 0; i < Mspikesr; i++)
	{
		vector<vector<double> > spikesd_bin;
		spikesd_bin.clear();
		int rb  = spikesr[i].first;
		for(int k = 0; k < Lfft; k++)
		{
			int km1 = k == 0 ? Lfft-1 : k-1;
			int kp1 = k == Lfft-1? 0 : k+1;

			if( YPM[rb][k].real() > YPM[rb][kp1].real() &&
					YPM[rb][k].real() > YPM[rb][km1].real() &&
					YPM[rb][k].real() > (double)Nrange/2.0 )
			{
				vector<double> spikesd;
				vector<double> obj;
				spikesd.clear();
				obj.clear();
//				cout << k << " " <<YPM[rb][km1].real()<<" "<<YPM[rb][k].real()
//						<< " " << YPM[rb][kp1].real() << endl;
				spikesd.push_back(k);
				spikesd.push_back(YPM[rb][km1].real());
				spikesd.push_back(YPM[rb][k].real());
				spikesd.push_back(YPM[rb][kp1].real());
				double amp, del_k;
				double z1 = sqrt(spikesd[1]);
				double z2 = sqrt(spikesd[2]);
				double z3 = sqrt(spikesd[3]);
				peakinterp(amp, del_k, z1, z2, z3);
				obj.push_back(amp);
				obj.push_back(rangep[rb]);

				double vel = ((spikesd[0] + del_k - 1)/Lfft -0.5)*vua;
				obj.push_back(vel);

				spikesd_bin.push_back(spikesd);
				targets.push_back(obj);
			}
		}
	}

	cout << "detected targets are:" << endl;
	int Mtargets = targets.size();
	int Ntargets = targets[0].size();
	for(int i = 0; i < Mtargets; i++)
	{
		for(int j = 0 ; j < Ntargets; j++)
		{
			cout <<" " << targets[i][j];
		}
		cout << endl;
	}


}

int main()
{
	vector< valarray<complex<double> > > y;
	radar a_radar;
	get_a_radar(a_radar);
	get_radar_signal(a_radar, y);
	process_radar_signal(a_radar, y);
}





#endif /* SIMULATOR_CPP_ */
