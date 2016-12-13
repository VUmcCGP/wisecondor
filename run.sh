#!/usr/bin/env bash
#set -e

WORKDIR=./

DATADIR=dataFiles
REFSAMPLEDIR=refSamples
TESTSAMPLEDIR=testSamples

PYTHON=python2.7
BINSIZE=50000
BINSIZEREF=250000

echo '
Please note:'
echo '	Copyright (C) 2016 VU University Medical Center Amsterdam
	Author: Roy Straver (github.com/rstraver)

	This file is part of WISECONDOR
	WISECONDOR is distributed under the following license:
	Attribution-NonCommercial-ShareAlike, CC BY-NC-SA
	(https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode)
	This license is governed by Dutch law and this license is subject
	to the exclusive jurisdiction of the courts of the Netherlands.
	'


# Turn bam files into npzs
function prepSamples {
	DIR=$1
	echo "	Preparing samples in $DIR"


	if [ `ls -1 $DIR/*.bam 2>/dev/null | wc -l` != 0 ];
	then
		for SAMPLE in `ls $DIR/*.bam`
		do
			if [ ! -f ${SAMPLE//.bam/.npz} ];
			then
				echo "		Working on: $SAMPLE -> ${SAMPLE//.bam/.npz}"
				$PYTHON wisecondor.py convert $SAMPLE ${SAMPLE//.bam/.npz} -binsize $BINSIZE > ${SAMPLE//.bam/_convert.log}
			fi
		done
	fi
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
	if [ ! -f $WORKDIR/$DATADIR/reference.npz ];
	then
		echo "	Missing: $DATADIR/reference"
		if [ `ls $WORKDIR/$REFSAMPLEDIR/* | wc -l` -eq 0 ];
		then
			echo ""
			echo "Please provide several *.bam files in $REFSAMPLEDIR/ obtained from healthy pregnancies. More samples yield better results."
			exit
		else
			prepSamples $WORKDIR/$REFSAMPLEDIR
			echo "		Found" `ls $WORKDIR/$REFSAMPLEDIR/*.npz | wc -l` "npz files in $REFSAMPLEDIR/"
			$PYTHON wisecondor.py newref $WORKDIR/$REFSAMPLEDIR/*.npz $WORKDIR/$DATADIR/reference.npz -binsize $BINSIZEREF
		fi
	else
		echo "	$DATADIR/reference found"
		echo "		Note: to update your reference with more samples remove $DATADIR/reference.npz"
	fi
	echo "		Note: previously converted .bam files can be removed from $REFSAMPLEDIR/ to save disk space"


# Testing samples itself
echo ""
echo "Testing data"
	if [ `ls $WORKDIR/$TESTSAMPLEDIR/* | wc -l` -eq 0 ];
	then
		echo ""
		echo "Please provide *.bam files in $TESTSAMPLEDIR/ to test pregnancies for aberrations."
		exit
	fi
	prepSamples $WORKDIR/$TESTSAMPLEDIR

	echo "	Testing samples"
	for SAMPLE in `ls $WORKDIR/$TESTSAMPLEDIR/*.npz | grep -v '_out.npz'`
	do
		if [ ! -f ${SAMPLE//.npz/_out.npz} ];
		then
			echo "		Working on: $SAMPLE"
			$PYTHON wisecondor.py test $WORKDIR/$SAMPLE $WORKDIR/${SAMPLE//.npz/_out.npz} $WORKDIR/$DATADIR/reference.npz > $WORKDIR/${SAMPLE//.npz/_test.log}
			$PYTHON wisecondor.py plot $WORKDIR/${SAMPLE//.npz/_out.npz} $WORKDIR/${SAMPLE//.npz/_plot} > $WORKDIR/${SAMPLE//.npz/_plot.log}
			$PYTHON wisecondor.py report $WORKDIR/$SAMPLE $WORKDIR/${SAMPLE//.npz/_out.npz} > $WORKDIR/${SAMPLE//.npz/.txt}
		fi
	done


echo ""
echo "Script finished"
