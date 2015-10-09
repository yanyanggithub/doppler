set terminal png size 1000, 500 enhanced font "Times 36"
set grid

set output "chirp_signal.png"
set title "Chirp Signal"
set xlabel "time (usec)"
set ylabel "amplitude (dB)"
plot "./chirp.dat" using 1:2 title 'real' w lines, \
	 "./chirp.dat" using 1:3 title 'imag' w lines
	 
reset	 
set output "y.png"
set key off
set title "Simulated Radar Signal (20 pulses)"
set xlabel "distance (km)"
set ylabel "amplitude (dB)"
plot for[i = 2:21] "./y.dat" using 1:(column(i)) w lines

reset
set output "ydb.png"
set key off
#set contour base
#set dgrid3d 40,40
unset surf
#set style line 1 lt 10 lw .5
set pm3d at s hidden3d 1
set title "FAST-TIME/SLOW-TIME PLOT OF RAW DATA"
set xlabel "pulse number"
set ylabel "range (km)"
splot "./ydb.dat" using 1:2:3 w lines

reset
set output "doppler_range.png"
set key off
set contour base
set cntrparam level incremental 46.6286, -5, 46.6286-30
set palette rgbformulae 33,13,10
unset surface
set style line 1 lt 6 lw 1.5
set pm3d at s hidden3d 1
set title "RANGE-DOPPLER PLOT OF UNPROCESSED DATA"
set xlabel "range"
set ylabel "velocity (km)"
splot "./doppler_range.dat" using 1:2:3 w lines

reset
set output "doppler_range_compressed.png"
set key off
set contour base
set cntrparam level incremental 72.494, -5, 72.494-20
set palette rgbformulae 33,13,10
unset surface
set style line 1 lt 6 lw 1.5
set pm3d at s hidden3d 1
set title "RANGE-DOPPLER PLOT OF PULSE-COMPRESSED DATA"
set xlabel "range"
set ylabel "velocity (km)"
splot "./doppler_range_compressed.dat" using 1:2:3 w lines

reset
set output "doppler_range_canceller.png"
set key off
set contour base
set cntrparam level incremental 81.6793, -5, 81.6793-25
set palette rgbformulae 33,13,10
unset surface
set style line 1 lt 6 lw 1.5
set pm3d at s hidden3d 1
set title "RANGE-DOPPLER PLOT OF PULSE-COMPRESSED DATA"
set xlabel "range"
set ylabel "velocity (km)"
splot "./doppler_range_canceller.dat" using 1:2:3 w lines



	
	 