#!/bin/bash
#

dir=$1
outdir=/media/ubuntu/data/dsa110/T1/${dir}/
python=/home/ubuntu/anaconda3/envs/casa/bin/python
filplot=/home/ubuntu/dana/code/dsa110-T3/dsaT3/filplot.py
mkdir -p $outdir

if ls /data/dsa110/T1/corr17/${dir}/fil_* 1> /dev/null 2>&1; then
    for fl in `ls -td /data/dsa110/T1/corr17/${dir}/fil_* | grep -v "done"`; do
	outname=$outdir$(basename $fl).done
	echo $outname
	if [ ! -f $outname ]; then
	    snum=`echo $fl | sed 's/\_/ /' | awk '{print $2}'`		
	    if ls /data/dsa110/T3/corr01/${dir}/*${snum}*.json 1> /dev/null 2>&1; then
		$python $filplot $dir $snum -s -c
		touch ${outname}
		echo $fl
	    fi
	fi
    done
fi
