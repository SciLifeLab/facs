1->Usage
a->'good_build' is for building bloom filter for reference file. It Supports large genome file (>4GB),for instance, human genome.
b->'simple_check' is for checking, also supporting large bloom file.
c->'simple_remove' is for removing.

2->How to Compile
gcc -c good_build.c -O3 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE -O3
gcc -c simple_check_1_ge.c -O3 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE -O3 -fopenmp
gcc -c simple_remove.c -O3 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE -O3 -fopenmp

gcc -o good_build good_build.o bloom.o suggestions.o lookup8.o -lm -O3 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE 
gcc -o simple_check simple_check_1_ge.o bloom.o suggestions.o lookup8.o -lm -O3 -fopenmp -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE
gcc -o simple_check simple_remove.o bloom.o suggestions.o lookup8.o -lm -O3 -fopenmp -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE 

3->Mannual 
After compiling, simply execute the scripts like './good_build' or './simple_remove' without any arguments will trigger an instruction page.
The whole query will be looked like these 
./good_build -m 1 -k 21 -e 0.005 -l tzcoolman2 (-p /home/user/desktop/)
./simple_check -m 1 -q test.fna -l tzcoolman3 -t 0.8 -s 1 (-p /home/user/desktop/)
./simple_remove -m 1 -q test.fna -l tzcoolman4 -t 1 (-p /home/user/desktop/)

4->Be advised
a-> All three scripts can be executed on both Linux and Mac system. But they dont support large bloom filter building and loading on MAC system.
b-> If you are going to run them on any clusters in Uppmax, e.g. Kalkyl, tintin... Make sure use that when you are building a filter for a large 
genome, or when you are loading a large filter, dont directly execute them, use sbatch script instead. 
Info can be found here http://www.uppmax.uu.se/support/user-guides/kalkyl-user-guide
c-> These scripts support fasta and fastq format. Make sure you use the correct extension name. 
For fasta files, please use .fna or .fasta.
For fastq files, please use .fastq.