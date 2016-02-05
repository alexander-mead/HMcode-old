reset

#This gnuplot script plots the output file 'power.dat' that is spat out of HMcode.
#Simply load up gnuplot (type gnuplot in the terminal) and then type "gnuplot>load 'plot.p'"
#The plot should then be the non-linear spectrum at 16 redshifts

set log x
set xrange [1e-3:1e4]
set xlabel 'k/(h Mpc^{-1})'

#set log y
#set yrange [1e-8:1e4]
#set ylabel '{/Symbol D}^2(k)'
#set format y '10^{%T}'
set yrange [0.98:1.02]

file1='power_old.dat'
file2='power_new.dat'

set key top left

col(i)=sprintf("#%1x%1x0000",i-1,i-1)

set title 'Test showing ratio between old and new versions of HMcode'

unset key

plot \
     0.99 w l lw 2 lt 2 lc rgb 'black',\
     1.01 w l lw 2 lt 2 lc rgb 'black',\
     for[i=1:16] '<paste '.file1.' '.file2.'' u 1:(column(i+1)/column(i+18)) w l lw 2 lc rgb col(i) noti



