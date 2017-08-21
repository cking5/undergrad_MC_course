awk -v param=$1: '{if($1=="timestamp:")printf("%f\t",$2);else if($1==param)printf("%f\n",$2);}' YAMLDATA.00* > $1.dat
