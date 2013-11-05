import os
import sys
import warnings

import ConfigParser

config = ConfigParser.SafeConfigParser()
conf_file = config.read([os.path.expanduser('~/.facsrc'), '.facsrc',
            'facs.conf', 'facs.cfg', '/etc/facs.conf'])
try:
# First config file found wins
    config.readfp(open(conf_file[0]))

    SERVER = config.get('facs', 'SERVER').rstrip()
    FACS_DB = config.get('facs', 'FACS_DB').rstrip()
    FASTQ_SCREEN_DB = config.get('facs', 'FASTQ_SCREEN_DB').rstrip()
    DECONSEQ_DB = config.get('facs', 'DECONSEQ_DB').rstrip()
    USERNAME = config.get('facs', 'USERNAME').rstrip()
    PASSWORD = config.get('facs', 'PASSWORD').rstrip()
    if config.has_option('facs', 'WAKE'):
        WAKE = bool(config.get('facs', 'WAKE'))
    else:
        WAKE = False 
except:
	warnings.warn("Please make sure you've created your own configuration file (i.e: ~/.facsrc) as stated in README.md")
