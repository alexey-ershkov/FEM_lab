set terminal png size 1100,750

set xlabel "x"
set ylabel "Error"

set output '../graphs/5_error.png'
plot '../cmake-build-debug/data/linear_error5' u 1:2 w l title 'Linear 5', \
'../cmake-build-debug/data/cubic_error5' u 1:2 w l title 'Cubic 5', \

set output '../graphs/20_error.png'
plot '../cmake-build-debug/data/linear_error20' u 1:2 w l title 'Linear 20', \
'../cmake-build-debug/data/cubic_error20' u 1:2 w l title 'Cubic 20', \

set output '../graphs/40_error.png'
plot '../cmake-build-debug/data/linear_error40' u 1:2 w l title 'Linear 40', \
'../cmake-build-debug/data/cubic_error40' u 1:2 w l title 'Cubic 40', \