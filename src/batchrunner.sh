#!/bin/bash

#
# Batchrunner. Puts a directory of inputs through the rr machine.
#

echo "batchrunner"

echo "Runs ridgerunner on a directory of inputs in an organized way."
echo "Checking for directories..."

if [ -d inputs ]; then

    echo "Found ./inputs directory."

else 

    echo "Need the directory ./inputs to exist."
    exit 1;

fi

if [ -d working ]; then

    echo "Found ./working directory."
    echo "Cleaning ./working..."

    cd ./working
    rm -fr *
    cd ..

else 

    echo "Creating ./working directory."
    mkdir ./working

fi

if [ -d trouble ]; then

    echo "Found ./trouble directory. (./trouble not cleaned)"
 
else 

    echo "Creating ./trouble directory."
    mkdir ./trouble

fi

if [ -d done ]; then

    echo "Found ./done directory. (./done not cleaned)"
 
else 

    echo "Creating ./done directory."
    mkdir ./done

fi

echo "Beginning runs..."

cd ./inputs

for basefile in *.vect
do
  
  basename=`echo $basefile | sed s/.vect//g`
  
  res=5

  thisbase=$basefile.r$res
  thisvect=$thisbase.vect
      
  echo cp $basefile ../working/$thisvect
  cp $basefile ../working/$thisvect
  cd ../working
      
  while [ $res -le 10 ]; do 

      rdrop=-1
      count=1

      echo "Running " $basename " at resolution " $res 
      
      oldvect=$thisvect
      thisbase=$basefile.r$res
      thisvect=$thisbase.vect

      cp $oldvect $thisvect

      while [ $rdrop -le -0.1 -a $count < 10 ]; do
      
	  echo "Running rr on " $thisvect " for 20,000 steps..."
	  ridgerunner -a -r $res $thisvect -s 20000 > $thisvect.out
	  echo "Run completed."

	  finalvect=./$thisbase.rr/$thisbase.final.vect

	  echo "Checking for " $finalvect  "..."

	  if [ -f $finalvect ]; then

	      echo "Found it."

	  else 

	      echo "The run seems to have failed. Moving to ./trouble..."
	      mv $thisbase.* ../trouble
	      break 2                     # break out of the drop and res loops

	  fi

	  echo "Checking to see whether ropelength dropped... "
	  startrop=`ropelength -q $thisvect`
          endrop=`ropelength -q $finalvect`
	  echo "     Initial ropelength " $startrop 
	  echo "     Final   ropelength " $endrop
	  echo "     Ropelength change  " $endrop - $startrop
	  rdrop=$endrop - $startrop

	  echo "Copying " $finalvect " to " $thisvect " for next run."
	  cp $finalvect ./$thisvect
	  count=$count+1
	  echo ""

      done

      # We are now ready to change resolution. We kick this copy into ./done

      echo "Finished " $thisvect " at resolution " $res 
      echo "Moving " `ls $thisbase*` " to ./done" 
      cp $thisvect ../done/$thisvect
      mv $thisbase* ../done/    
	
  done

  cd ../inputs

done
