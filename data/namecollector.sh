#!/bin/bash
echo "datalist.txt. Create this file with ./namecollector > datalist.txt"
N="0"
for name in $( ls *.vect ); do 
   N=$(( $N + 1 )); 
   NARRAY[$N]=$name;
   if [ $N -eq 5 ] 
       then 

       echo ${NARRAY[0]} ${NARRAY[1]} ${NARRAY[2]} ${NARRAY[3]} ${NARRAY[4]} ${NARRAY[5]} " \\"
       N="0"
   fi
done
