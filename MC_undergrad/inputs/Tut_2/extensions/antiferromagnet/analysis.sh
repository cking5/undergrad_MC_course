#!/bin/bash
# r.j.grant@bath.ac.uk
# Script to strip Magnetisation from DLMONTE 64x64 2D ISING PTFILE and produce:
# 1. a time sequence of magnetisation
# 2. a histogram of data
# 3. calculate the magnetisation
# 4 plot 1 and 2
awk '/system temperature/{print $5}' OUTPUT.000 > T.temp
grep --text -A14 'averages and fluctuations' OUTPUT.000 | tail -1 > output.temp
awk '{print $2}' output.temp > output2.temp
awk '{print}' T.temp output2.temp > output
sed 'N;s/\n/ /' output > energy.dat
awk '{if((NR + 5)%11==6) phi=$1};{if(NR%11==6)print phi, ($1-$2);}' PTFILE.000 > M_seq.dat
awk '{hist[$2+4097]++;} END{for(i=1;i<8194;i++)print (i-4097)/4096,hist[i]/NR;}' M_seq.dat > M_hist.dat
awk '{M+=$2*$2;} END {print sqrt(M/NR)/4096;}' M_seq.dat > Mm.dat
awk '{print}' T.temp Mm.dat > Mn.dat
sed 'N;s/\n/ /' Mn.dat > M.dat
rm T.temp output.temp output2.temp output Mm.dat Mn.dat
gnuplot plot.gpl
