"""Compare performance against fastq_screen.

Given previous test results stored in an external database (specified in your
~/.facsrc file), downloads this results and plot them.
"""
import logbook
import sys
import couchdb

from facs.utils import config

log = logbook.Logger('FACS')

def _fetch_results(couch, database):
    """Fetch all documents a from database
    """
    log.info("Fetching all documents from database %s" % database)
    docs = []
    [docs.append(couch[database].get(doc)) for doc in couch[database]]
    log.info("Fetched %s documents" % str(len(docs)))
    return docs


def facs_vs_fastq_screen():
    stream = logbook.StreamHandler(sys.stdout, level=logbook.INFO)
    with stream.applicationbound():
        log.info("Establishing connection with database %s" % config.SERVER)
        couch = couchdb.Server(config.SERVER)
        couch.resource.credentials = (config.USERNAME, config.PASSWORD)
        facs_results = _fetch_results(couch, config.FACS_DB)
        fastq_screen_results = _fetch_results(couch, config.FASTQ_SCREEN_DB)
