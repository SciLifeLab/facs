#!/bin/bash
#SBATCH -p node -n 8
#SBATCH -t 02:00:00
#SBATCH -J facs
#SBATCH -o log.out
#SBATCH -e log.err
#SBATCH -D /bubo/home/h5/roman/dev/facs/tests
#SBATCH -A a2010002
#SBATCH --qos=seqver
#SBATCH -C thin
#SBATCH --mail-type=all
#SBATCH --mail-user=roman@scilifelab.se

nosetests --with-timer -v -s test_innocentive.py
mv log.out /proj/b2012094/Innocentive/data/logs/`date +%Y-%m-%d`.log
mv log.err /proj/b2012094/Innocentive/data/logs/`date +%Y-%m-%d`.err
