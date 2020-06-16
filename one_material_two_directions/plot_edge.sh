source activate mp1.8

interp=49
echo "$interp" 
# mpirun -np $nworkers python3 band.py $interp > band.out
python3 band.py $interp > band.out
echo "Calculation finished"
grep tefreqs: band.out > bandte.dat
# grep tmfreqs: band.out > bandtm.dat
python3 plot_edge_band.py $interp
