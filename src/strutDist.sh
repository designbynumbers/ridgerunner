#!/bin/sh

cp /tmp/struts.vect /tmp/foo
lines=`head -n 2 /tmp/foo | tail -n 1 | awk '{print $1}'`
echo $lines struts
tail -n $lines /tmp/foo | sed "/,0.000000,0.000000,0/s///" | sort -n > /tmp/dat
echo "set key off" > /tmp/cmd
echo "plot \"/tmp/dat\" w lines, \"/tmp/dat\" w points" >> /tmp/cmd
echo "q" >> /tmp/cmd
gnuplot < /tmp/cmd
cat /tmp/dat
rm -f /tmp/cmd
rm -f /tmp/dat
