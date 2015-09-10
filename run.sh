##############################################################################
#                                                                            #
#    Test a maternal DNA sample for fetal Copy Number Aberrations.           #
#    Copyright(C) 2013  TU Delft & VU University Medical Center Amsterdam    #
#    Author: Roy Straver, r.straver@vumc.nl                                  #
#                                                                            #
#    This file is part of WISECONDOR.                                        #
#                                                                            #
#    WISECONDOR is free software: you can redistribute it and/or modify      #
#    it under the terms of the GNU General Public License as published by    #
#    the Free Software Foundation, either version 3 of the License, or       #
#    (at your option) any later version.                                     #
#                                                                            #
#    WISECONDOR is distributed in the hope that it will be useful,           #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of          #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
#    GNU General Public License for more details.                            #
#                                                                            #
#    You should have received a copy of the GNU General Public License       #
#    along with WISECONDOR.  If not, see <http://www.gnu.org/licenses/>.     #
#                                                                            #
##############################################################################



set -e

WORKDIR=./

DATADIR=dataFiles
REFSAMPLEDIR=refSamples
TESTSAMPLEDIR=testSamples


# Obtain required packages
PACKAGELIST="python2.7 python-numpy python-matplotlib python-biopython samtools" #build-essential zlib1g-dev libncurses5-dev libncursesw5-dev" # default-jre"
INSTALLEDCOUNT=`dpkg -s $PACKAGELIST 2>/dev/null | grep 'ok installed' | wc -l`
PACKAGECOUNT=`echo $PACKAGELIST | wc -w`

echo ""
echo 'Please note this script relies on the following packages:' 
echo $PACKAGELIST
echo 'These are copyrighted by their respective owners.'
echo ""

if  [ ! $INSTALLEDCOUNT -eq $PACKAGECOUNT ]
then
	echo 'Additional packages are required'
	sudo apt-get install $PACKAGELIST
else
	echo 'All required packages were found'
fi



# Turn bam files into pickles, pickles into gcc
function prepSamples {
	DIR=$1
	echo "	Preparing samples in $DIR"
	
	for SAMPLE in `ls $DIR/*.bam`
	do
		if [ ! -f ${SAMPLE//.bam/.pickle} ]; 
		then
			echo "		Working on: $SAMPLE -> ${SAMPLE//.bam/.pickle}"
			samtools rmdup -s $SAMPLE - | samtools view - -q 1 | python2.7 consam.py -outfile ${SAMPLE//.bam/.pickle}
		fi
	done
	
	for SAMPLE in `ls $DIR/*.pickle`
	do
		if [ ! -f ${SAMPLE//.pickle/.gcc} ]; 
		then
			echo "		Working on: $SAMPLE -> ${SAMPLE//.pickle/.gcc}"
			python2.7 ./gcc.py $SAMPLE ./$DATADIR/gccount ${SAMPLE//.pickle/.gcc} > ${SAMPLE//.pickle/.log}
		fi
	done
}



# Create working directories
echo ""
echo "Directory structure"
for DIRECTORY in $DATADIR $REFSAMPLEDIR $TESTSAMPLEDIR
do
	if [ ! -d "$WORKDIR/$DIRECTORY" ]; 
	then
		echo "	Creating directory: $DIRECTORY"
		mkdir $WORKDIR/$DIRECTORY
	else
		echo "	$DIRECTORY found"
	fi
done



# Check for reference files
echo ""
echo "Reference data"
	if [ ! -f $WORKDIR/$DATADIR/gccount ]; 
	then
		echo "	Missing: $DATADIR/gccount"
		
		if [ ! -f $WORKDIR/$DATADIR/*.fa* ]; 
		then
			echo ""
			echo "Please provide a single *.fa* file in $DATADIR/ containing the reference used for mapping reads." 
			exit
		else
			echo "Creating gccount"
			python2.7 countgc.py $WORKDIR/$DATADIR/*.fa* $WORKDIR/$DATADIR/gccount
		fi

	else
		echo "	$DATADIR/gccount found"
		if [ -f $WORKDIR/$DATADIR/*.fa ]; 
		then
			echo "		Note: you may (re)move $DATADIR/*.fa* to save disk space"
		fi
	fi
	
	if [ ! -f $WORKDIR/$DATADIR/reference ]; 
	then
		echo "	Missing: $DATADIR/reference"
		if [ `ls $WORKDIR/$REFSAMPLEDIR/*.bam | wc -l` -eq 0 ]; 
		then
			echo ""
			echo "Please provide several *.bam files in $REFSAMPLEDIR/ obtained from healthy pregnancies. More samples yield better results."		
			exit	
		else
			prepSamples $WORKDIR/$REFSAMPLEDIR
			echo "		Found" `ls $WORKDIR/$REFSAMPLEDIR/*.pickle | wc -l` "pickle files in $REFSAMPLEDIR/"
			python2.7 newref.py  $WORKDIR/$REFSAMPLEDIR $WORKDIR/$DATADIR/reference
		fi
	else
		echo "	$DATADIR/reference found"
		echo "		Note: to update your reference with more samples remove $DATADIR/reference"
	fi
	echo "		Note: previously converted .bam files can be removed from $REFSAMPLEDIR/ to save disk space"
	


# Testing samples itself
echo ""
echo "Testing data"
	if [ `ls $WORKDIR/$TESTSAMPLEDIR/*.bam | wc -l` -eq 0 ]; 
	then
		echo ""
		echo "Please provide *.bam files in $TESTSAMPLEDIR/ to test pregnancies for aberrations."	
		exit
	fi
	prepSamples $WORKDIR/$TESTSAMPLEDIR
	
	echo "	Testing samples"
	for SAMPLE in `ls $WORKDIR/$TESTSAMPLEDIR/*.gcc`
	do
		if [ ! -f ${SAMPLE//.gcc/.plot} ]; 
		then
			echo "		Working on: $SAMPLE"
			python2.7 test.py $SAMPLE $WORKDIR/$DATADIR/reference $WORKDIR/${SAMPLE//.gcc/.plot} > ${SAMPLE//.gcc/.txt}
			python2.7 plot.py ${SAMPLE//.gcc/.plot} ${SAMPLE//.gcc/} >> ${SAMPLE//.gcc/.txt}
			python2.7 makejson.py ${SAMPLE//.gcc/.plot} > ${SAMPLE//.gcc/.json}
		fi
	done
	
	
	
echo ""
echo "Script finished" 
