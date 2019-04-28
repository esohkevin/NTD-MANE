
set terminal png nocrop enhanced size 640,480 font "arial,12.0"
set output 'poten_energy.png'
set format y '%5.1fe'
set xlabel "Energy Minimization Step"
set ylabel "Potential Energy (kJ mol-1)"
plot "potential.xvg" using 1:2 title 'Potential' with lines

