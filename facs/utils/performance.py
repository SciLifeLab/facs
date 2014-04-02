"""Compare performance against fastq_screen.

Given previous test results stored in an external database (specified in your
~/.facsrc file), downloads this results and plot them.

XXX: Needs serious refactoring, way too much code repetition :/
"""
import logbook
import os
import glob
import sys
import json
import couchdb
import datetime
from collections import defaultdict
from itertools import izip

from facs.utils import config

log = logbook.Logger('FACS')

def _fetch_results(couch, database):
    """Fetch all documents a from database
    """
    log.info("Fetching all documents from database %s" % database)
    docs = []
    db = couch[database]
    [docs.append(db.get(doc)) for doc in db]
    log.info("Fetched %s documents" % str(len(docs)))
    return docs

def facs_vs_deconseq():
    """ Compare FACS against bwa-based deconseq
    """
    facs, deconseq = "", ""
    results = defaultdict(list)

    with open("facs.json") as fh:
        facs = json.load(fh)

    with open("deconseq.json") as fh:
        deconseq = json.load(fh)

    for decon in deconseq:
        if decon.get('sample'):
            for fcs in facs:
                if fcs.get('sample'):
                    if os.path.basename(fcs.get('sample')) == os.path.basename(decon.get('sample')):
                        if fcs.get('begin_timestamp') and decon.get('start_timestamp'):
                            begin_fcs = fcs['begin_timestamp']
                            end_fcs = fcs['end_timestamp']

                            begin_deco = decon['start_timestamp']
                            end_deco = decon['end_timestamp']

                            # remove the UTC offset (+0200) (%z does not parse it out)
                            # http://docs.python.org/2/library/datetime.html#strftime-strptime-behavior
                            begin_fcs = begin_fcs[:-5]
                            end_fcs = end_fcs[:-5]

                            dt_b_fcs = datetime.datetime.strptime( begin_fcs, "%Y-%m-%dT%H:%M:%S.%f" )
                            dt_e_fcs = datetime.datetime.strptime( end_fcs, "%Y-%m-%dT%H:%M:%S.%f" )

                            # How log did FACS run?
                            delta_fcs = dt_e_fcs - dt_b_fcs


                            # Remove the final 'Z' in timestamp
                            begin_deco = begin_deco[:-1]
                            end_deco = end_deco[:-1]

                            dt_b_deco = datetime.datetime.strptime( begin_deco, "%Y-%m-%d %H:%M:%S.%f" )
                            dt_e_deco = datetime.datetime.strptime( end_deco, "%Y-%m-%d %H:%M:%S.%f" )

                            # How long did deconseq run?
                            delta_deco = dt_e_deco - dt_b_deco

                            # Fetch contamination rates for both
                            contam_fcs = fcs.get('contamination_rate')
                            contam_deco = decon.get('contamination_rate')

                            results[run] = dict(delta = delta_deco.total_seconds(),
                                                contam_deco = contam_deco,
                                                delta_facs = delta_fcs.total_seconds(),
                                                contam_facs = contam_fcs,
                                                sample_deco = decon.get('sample'),
                                                sample_facs = fcs.get('sample'),
                                                filter_deco = decon.get('reference'),
                                                filter_facs = fcs.get('bloom_filter'))

    return json.dumps(results)

def facs_vs_fastq_screen():
    """ Work from the json files on disk instead of fetched from DB
    """
    facs, fastq_screen = "", ""
    results = defaultdict(list)

    with open("facs.json") as fh:
        facs = json.load(fh)

    with open("fastq_screen.json") as fh:
        fastq_screen = json.load(fh)

    for run, (fqscr, fcs) in enumerate(izip(fastq_screen, facs)):
        if fqscr.get('sample') and fcs.get('sample'):
            if os.path.basename(fcs['sample']) == fqscr['sample']:
                facs_filt_name, _ = os.path.splitext(os.path.basename(fcs['bloom_filter']))
                if facs_filt_name == fqscr['fastq_screen_index']:
                    # Do not assume the test run went well
                    if len(fqscr['organisms']) > 0 and fqscr.get('memory_usage', None) and \
                            fcs.get('memory_usage', None):
                        if fqscr.get('begin_timestamp') and fcs.get('begin_timestamp'):
                        # Fetch timing info for each program
                            begin = fqscr['begin_timestamp']
                            end = fqscr['end_timestamp']

                            dt_b = datetime.datetime.strptime( begin, "%Y-%m-%d %H:%M:%S.%fZ" )
                            dt_e = datetime.datetime.strptime( end, "%Y-%m-%d %H:%M:%S.%fZ" )

                            # Fastqscreen runtime
                            delta = dt_e - dt_b

                            begin_fcs = fcs['begin_timestamp']
                            end_fcs = fcs['end_timestamp']

                            # remove the UTC offset (+0200) (%z does not parse it out)
                            # http://docs.python.org/2/library/datetime.html#strftime-strptime-behavior
                            begin_fcs = begin_fcs[:-5]
                            end_fcs = end_fcs[:-5]

                            dt_b_f = datetime.datetime.strptime( begin_fcs, "%Y-%m-%dT%H:%M:%S.%f" )
                            dt_e_f = datetime.datetime.strptime( end_fcs, "%Y-%m-%dT%H:%M:%S.%f" )

                            # FACS runtime
                            delta_fcs = dt_e_f - dt_b_f

                            # Fetch contamination rates for both
                            contam_fcs = fcs.get('contamination_rate')
                            contam_fqscr = fqscr.get('contamination_rate')

                            # Memory usage
                            mem_facs = fcs.get('memory_usage')
                            mem_fqscr = fqscr.get('memory_usage')

                            results[run] = dict(delta = delta.total_seconds(),
                                                contam_fqscr = contam_fqscr,
                                                delta_facs = delta_fcs.total_seconds(),
                                                contam_facs = contam_fcs,
                                                sample_fqscr = fqscr.get('sample'),
                                                sample_facs = fcs.get('sample'),
                                                filter_fqscr = fqscr.get('fastq_screen_index'),
                                                filter_facs = fcs.get('bloom_filter'),
                                                threads_facs = fcs.get('threads'),
                                                threads_fqscr = fqscr.get('threads'),
                                                mem_facs = mem_facs,
                                                mem_fqscr = mem_fqscr)


    return json.dumps(results)


def fetch_couchdb_results():

    stream = logbook.StreamHandler(sys.stdout, level=logbook.INFO)
    with stream.applicationbound():
        log.info("Establishing connection with database %s" % config.SERVER)

        couch = couchdb.Server(config.SERVER)
        couch.resource.credentials = (config.USERNAME, config.PASSWORD)

        facs_results = _fetch_results(couch, config.FACS_DB)
        fastq_screen_results = _fetch_results(couch, config.FASTQ_SCREEN_DB)
        #deconseq_results = _fetch_results(couch, config.DECONSEQ_DB)

        with open("facs.json", 'w') as fh:
            json.dump(facs_results, fh)

        with open("fastq_screen.json", 'w') as fh:
            json.dump(fastq_screen_results, fh)

        #with open("deconseq.json", 'w') as fh:
        #    json.dump(deconseq_results, fh)

if __name__ == "__main__":

    # Fetch CouchDB dumps locally for now
    # iriscouch is not the fastest nor most reliable DB around
    if not glob.glob('*.json'):
        fetch_couchdb_results()

    log.info("Comparing runtimes of FACS vs fastq_screen")
    print facs_vs_fastq_screen()

    #log.info("Comparing runtimes of FACS vs deconseq")
    #facs_vs_deconseq()
