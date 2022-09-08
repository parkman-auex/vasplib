set encoding iso_8859_1
#set terminal  postscript enhanced color
#set output 'surfdos_r.eps'
#set palette rgbformulae 33,13,10
set style data linespoints
unset ztics
unset key
set pointsize 0.8
set pm3d
set border lw 3
set size 0.8, 1
set origin 0.1, 0
#set size square
#set view equal xyz
set view map
#set cbtics font ",48"
#set xtics font ",48"
#set ytics font ",48"
#set ylabel font ",48"
set ylabel "Energy (eV)"
#set xtics offset 0, -1
#set ylabel offset -6, 0 
set xrange [0:            0.86499]
set yrange [          -0.20000:           0.19960]
set xtics ("Z"  0.00000,"G"  0.43250,"Z"  0.86499)
set arrow from  0.43250,  -0.20000 to  0.43250,   0.19960 nohead front lw 1
set pm3d interpolate 2,2
unset cbtics
splot 'dos.dat_r' u 1:2:3 w pm3d
