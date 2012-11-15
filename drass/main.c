#include "check.h"
#include "build.h"
#include "remove.h"
#include "remove_l.h"
#include "big_query.h"
/*------------------------------*/  
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
/*------------------------------*/ 
int k_mer, help;
float tole_rate, error_rate, sampling_rate;
char *prefix, *list, *ref, *source, *mode;

/*------------------------------*/ 
//void descript ();
//void init (int argc, char **argv);
/*------------------------------*/ 

/*------------------------------*/ 
void
descript () {
  printf
    ("##########################################################################\n");
  printf ("#NAME\n");
  printf ("#DRASS--> Rapid and Accurate Sequence Decontamination System\n");
  printf ("# Contamination Removal tool\n");
  printf ("#\n");
  printf
    ("# extended from FACS (Fast and Accurate Classification of Sequences)\n");
  printf ("#\n");
  printf ("# facs.scilifelab.se\n");
  printf ("#\n");
  printf ("# Author Enze Liu\n");
  printf ("#\n");
  printf ("# enze.liu@scilifelab.se\n");
  printf
    ("##########################################################################\n");
  printf ("\n");
  printf ("DESCRIPTION\n");
  printf ("\n");
  printf
    ("    Bloom filters are a lightweight duplicate detection algorithm proposed\n");
  printf ("    by Burton Bloom\n");
  printf
    ("    (http://portal.acm.org/citation.cfm?id=362692&dl=ACM&coll=portal), with\n");
  printf
    ("    applications in stream data processing, among others. Bloom filters are\n");
  printf
    ("    a very cool thing. Where occasional false positives are acceptable,\n");
  printf
    ("    bloom filters give us the ability to detect duplicates in a fast and\n");
  printf ("    resource-friendly manner.\n");
  printf ("\n");
  printf
    ("    The allocation of memory for the bit vector is handled in the c layer,\n");
  printf
    ("    but perl's oo capability handles the garbage collection. when a\n");
  printf
    ("   Bloom::Faster object goes out of scope, the vector pointed to by the c\n");
  printf
    ("    structure will be free()d. to manually do this, the DESTROY builtin\n");
  printf ("    method can be called.\n");
  printf ("\n");
  printf
    ("    A bloom filter perl module is currently avaible on CPAN, but it is\n");
  printf
    ("    profoundly slow and cannot handle large vectors. This alternative uses a\n");
  printf
    ("    more efficient c library which can handle arbitrarily large vectors (up\n");
  printf ("    to the maximum size of a long long datatype (at least\n");
  printf ("    9223372036854775807, on supported systems ).\n");
  printf ("\n");
  printf
    ("    FACS is a novel algorithm, using Bloom filter to accurately and rapidly\n");
  printf
    ("    align sequences to a reference sequence. FACS was first optimized and\n");
  printf
    ("    validated using a synthetic metagenome dataset. An experimental metagenome\n");
  printf
    ("    dataset was then used to show that FACS is at least three times faster and\n");
  printf
    ("    more accurate than BLAT and SSAHA2 in classifying sequences when using\n");
  printf ("    references larger than 50Mbp.\n");
  printf ("\n");
}


void 
init (int argc, char **argv) 
{
    if (argc == 1 || !strcmp (argv[1], "-h") || !strcmp (argv[1], "-help")) {
        descript ();
        exit (0);
    }
     
    /*-------defaults-------*/ 
      
    k_mer = 21;
    tole_rate = 0.8;
    error_rate = 0.0005;
    sampling_rate = 1;

    help = 0;
      
    prefix = NULL;
    list = NULL;
    ref = NULL;
    source = NULL;
    mode = NULL;
      
    /*-------default-------*/ 
      
    int x;
      
    while ((x = getopt (argc, argv, "e:k:m:t:o:r:q:s:l:b:")) != -1) {

        switch (x) {
            case 'e':
            (optarg) && ((error_rate = atof (optarg)), 1);
            break;
            case 'k':
            (optarg) && ((k_mer = atoi (optarg)), 1);
            break;
            case 'm':
            (optarg) && ((mode = optarg), 1);
            break;
            case 't':
            (optarg) && ((tole_rate = atof (optarg)), 1); 
            break;
            case 's': 
            (optarg) && ((sampling_rate = atof (optarg)), 1);
            break;
            case 'o':	  
            (optarg) && ((prefix = optarg), 1);
            break;
            case 'r':
            (optarg) && ((ref = optarg), 1); 
            break;
            case 'q':  
            (optarg) && (source = optarg, 1);  
            break;
            case 'l':
            (optarg) && (list = optarg, 1);  
            break;
            case 'b':
            (optarg) && (help = atoi(optarg), 1);
            break;
            case '?':
            printf ("Unknown option: -%c\n", (char) optopt); 
            exit (0);
        } 
    } 

    if (help == 0)
       if ((source==NULL) && ((ref==NULL) && (list==NULL)))
        {
          perror ("Query or reference missing..."); 
          exit (0);
        }
}


int main (int argc, char **argv) 
{
  
init (argc, argv);
  
switch (*mode)
    {
	case 'b': 
	build_main (k_mer, error_rate, source, list, prefix,argv[0],help);
	break;
	case 'c':  
	check_main (source,ref,tole_rate,sampling_rate,list,prefix,help);  
	break;  
	case 'r':    
	remove_main (tole_rate, source, ref, list, prefix,help);   
	break; 
	/*case 'l':     
	remove_main_l (tole_rate, source, ref, list, prefix,help);     
	break;
    
    XXX
    main.o: In function `main':
    drass/main.c:182: undefined reference to `remove_main_l'
    */
    case 'g':
    bq_main (source,ref,tole_rate,sampling_rate,list,prefix,help);
    break;
	default:     
	printf ("Wrong mode: stupid idiot\n");    
	descript();
	break;
    }
  
return 0;

}
