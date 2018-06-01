reset

#set term aqua dashed

#This gnuplot script plots the output file 'power.dat' that is spat out of HMcode.
#Simply load up gnuplot (type gnuplot in the terminal) and then type "gnuplot>load 'plot.p'"
#The plot should then be the non-linear spectrum at 16 redshifts

#x-axis stuff
kmin=1e-3
kmax=1e2
set log x
set xrange [kmin:kmax]
set xlabel 'k/(h Mpc^{-1})'
set mxtics 10

#y-axis stuff
pmin=1e-8
pmax=1e4
set log y
set yrange [pmin:pmax]
set ylabel '{/Symbol D}^2(k)'
set format y '10^{%T}'
set mytics 1

#File to plot
file='power.dat'

#Colour box stuff
unset colorbox

#Key stuff
unset key

#Number of redshifts
n=16

#Now do the actual plotting
plot for[i=1:n] file u 1:(column(i+1)):(real(i-1)/real(n)) w l lw 2 dt 1 lc palette noti



