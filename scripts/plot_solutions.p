set terminal png size 1100,750

set output '../graphs/20.png'
set xlabel "x"
set ylabel "u" rotate by 0
set key bot right

set output '../graphs/5.png'
plot '../cmake-build-debug/data/analytic5' u 1:2 w l title 'Analytic 5', \
'../cmake-build-debug/data/linear5' u 1:2 w l title 'Linear 5', \
'../cmake-build-debug/data/cubic5' u 1:2 w l title 'Cubic 5', \

set output '../graphs/20.png'
plot '../cmake-build-debug/data/analytic20' u 1:2 w l title 'Analytic 20', \
'../cmake-build-debug/data/linear20' u 1:2 w l title 'Linear 20', \
'../cmake-build-debug/data/cubic20' u 1:2 w l title 'Cubic 20', \

set output '../graphs/40.png'
plot '../cmake-build-debug/data/analytic40' u 1:2 w l title 'Analytic 40', \
'../cmake-build-debug/data/linear40' u 1:2 w l title 'Linear 40', \
'../cmake-build-debug/data/cubic40' u 1:2 w l title 'Cubic 40', \