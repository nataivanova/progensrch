#!/usr/bin/env python3

import re
import numpy as np
import os
import pickle
import json
import logging

class ProgenSearchException(Exception):
    errors = []

    def __init__(self, errors):
        self.errors = errors

    def __str__(self):
        return str(self.errors)

class ProgenDBInitException (ProgenSearchException):
    pass

class QueryParametersException(ProgenSearchException):
    pass

class NoProgensFoundException(ProgenSearchException):
    pass


class ProgenitorQuery:

    # m1 -> Donor Mass {Msol}
    # m2 -> Accretor Mass {Msol}
    # mt -> log10(MT Rate) {Msol/yr}
    # p -> Orbital Period {days}
    # teff -> log10(Donor Effective Temp) {K}
    # labels = ['m1','m2','mt','p','teff']

    query = {}
    content = ""


    def __init__(self, qpath) -> None:

        logger = logging.getLogger('progentool')
        logger.info ('parsing the input file' + qpath)

        try:
            with open(qpath) as infile:
                self.content = infile.read()

                m = re.match (
                    r"""^(?P<isBH>[01])                                          \s*(\#.*)?\s*    # first line: 1 = BH
                      (?P<donor_m>    \d+\.?\d*(?: [ed][+-]?\d+)*)   \s*,\s* (?P<donor_M>    \d+\.?\d*(?: [ed][+-]?\d+)*) \s*(\#.*)?\s*    # donor mass
                      (?P<accretor_m> \d+\.?\d*(?: [ed][+-]?\d+)*)   \s*,\s* (?P<accretor_M> \d+\.?\d*(?: [ed][+-]?\d+)*) \s*(\#.*)?\s*    # accretor mass
                      (?P<mt_m>     -?\d+\.?\d*(?: [ed][+-]?\d+)*)   \s*,\s* (?P<mt_M>     -?\d+\.?\d*(?: [ed][+-]?\d+)*) \s*(\#.*?)\s*    # MT rate
                      (?P<period_m>   \d+\.?\d*(?: [ed][+-]?\d+)*)   \s*,\s* (?P<period_M>   \d+\.?\d*(?: [ed][+-]?\d+)*) \s*(\#.*?)\s*    # orbital period
                      (?P<teff_m>     \d+\.?\d*(?: [ed][+-]?\d+)*)   \s*,\s* (?P<teff_M>     \d+\.?\d*(?: [ed][+-]?\d+)*) \s*(\#.*)?\s*    # effective T
                      .*"""
                    , self.content
                    , flags = re.X)

        except Exception as e:
            logger.error("Error reading input file")
            logger.error(str(e))

        if m:
            self.query = { 'bhns': int(m.group('isBH'))
                           , 'm1': [float(m.group('donor_m')),    float(m.group('donor_M'))]
                           , 'm2': [float(m.group('accretor_m')), float(m.group('accretor_M'))]
                           , 'mt': [float(m.group('mt_m')),       float(m.group('mt_M'))]
                           , 'p':  [float(m.group('period_m')),   float(m.group('period_M'))]
                           , 'teff': [float(m.group('teff_m')),   float(m.group('teff_M'))]
                          }
        else:
            raise QueryParametersException('input file does not match template\n' + self.content)

        self.validate()

        return None

    def __str__(self) -> None:
        return str(self.query)

    def validate(self):
        """ Rudimentary input checks """
        query = self.query
        errors = []
        if (query['m1'][0] < 0.0 or query['m1'][0] > query['m1'][1]):
            errors.append("Wrong donor mass range")
        if (query['m2'][0] < 0.0 or query['m2'][0] > query['m2'][1]):
            errors.append("Wrong BH mass range")
        if (query['mt'][0] < -50.0 or query['mt'][0] > query['mt'][1]):
            errors.append("Wrong MT rate")
        if (query['p'][0] < 0.0 or query['p'][0] > query['p'][1]):
            errors.append("Wrong orbital period")
        if query['teff'][0] > query['teff'][1]:
            errors.append("Wrong Donor Teff")
        if (errors):
            raise QueryParametersException(errors)



class ProgenitorDatabase:
    """ Instanciate this class to initialize the interface to the database.
        In the initial implementation, the database is a collection of flat files.
        The constructor reads in all avaliable paths and infers the start/end masses
        for donor accretor for each file. Since donors only lose and accreators only
        gain mass, the mapping file path -> [mass ranges] works as a pair of indices
        on the two mass parameters.
    """

    logger = logging.getLogger ('progentool')

    db = {}

    def __init__(self, db_location: str) -> None:
        """ Read the file system structure and set up the mapping
            file name -> mass ranges
            A cache is set up and pickled to the file system,
            and is used to initialize the in-memory db by default.
        """

        self.logger.info("initializing the database")

        cache_location = db_location + '/db.pickle'
        try:
            with open(cache_location, 'rb') as cached_db:
                self.db = pickle.load(cached_db)
                self.logger.warning ('loaded cached database. If this is not desired, delete the file'
                                     + cache_location + ' to forcefully reinitialize the cache.')
        except Exception as e:
            self.logger.warning ('no db cache found. Reading db files from disk and initializing the cache. '
                                + 'This may take a while')
            self.init_cache(db_location)
            self.save_cache(cache_location)

        self.logger.info("loaded database with " + str(len(self.db)) + " files")


    def save_cache(self, cache_location: str) -> None:
        try:
            with open(cache_location, 'wb') as cached_db:
                pickle.dump (self.db, cached_db)
        except Exception as e:
            self.logger.error('cannot write the cache to location '
                              + cache_location
                              + "\n Will endeavour to proceed")

    def init_cache(self, db_location: str) -> None:
        try:
            for root, dirs, files in os.walk(db_location, followlinks=True):
                self.logger.debug ('walking the file system. We are at ' + root
                              + 'and see dirs '    + str(dirs)
                              + 'as well as files' + str(files) )

                for file in files:
                    fullpath = root + '/' + file
                    self.logger.debug('pondering over file ' + fullpath)

                    f = re.match(  r"""m_(?P<m1>\d+\.?\d+)_p_(?P<p>-?\d+\.?\d+)\.data$"""
                                   , file
                                   , flags = re.X)
                    g = re.match(  r"""
                    /.*
                    / (?P<isbh>[01])?        # 0 = NS, 1 = BH
                    /(?P<m2>[^/]*)?          # for black holes, the directory structure is split by accretor mass
                    /[ls]m[ls]p"""
                                   , root
                                   , flags = re.X)

                    if f:
                        if not g:
                            raise ProgenDBInitException('found a data file but unable to determine the system type')

                        self.logger.debug('adding file ' + file + ' to the database')
                        self.db[fullpath] = { 'isbh':  int(g.group('isbh'))
                                            , 'm1' : float(f.group('m1'))
                                            , 'p'  : float(f.group('p')) }
                        if self.db[fullpath]['isbh']:
                            self.db[fullpath]['m2'] = float(g.group('m2'))


        except Exception:
            raise ProgenDBInitException('cannot initialize the database')


        self.logger.info('first start through building the database, learned about ' + str(len(self.db)) + " files.\n Now caching some data..." )

        nomt = []

        for dbfile in self.db.keys():
            self.logger.debug ('scanning file ' + dbfile + ' to extract values')
            with open(dbfile, 'rb') as infile:
                donor_masses, accretor_masses, mt_rates, periods, teffs, ages, radii, dts  = pickle.load(infile)

            # Donors only lose mass, so:
            self.db[dbfile]['m1_min'] = donor_masses[-1]
            self.db[dbfile]['m1_max'] = donor_masses[0]

            # Conversely, accretors only gain mass:
            self.db[dbfile]['m2_min'] = accretor_masses[0]
            self.db[dbfile]['m2_max'] = accretor_masses[-1]

            # MT min/max:
            self.db[dbfile]['mt_min'] = min(mt_rates)
            self.db[dbfile]['mt_max'] = max(mt_rates)

            # period min/max:
            self.db[dbfile]['p_min'] = min(periods)
            self.db[dbfile]['p_max'] = max(periods)

            # Teff min/max
            self.db[dbfile]['t_min'] = min(teffs)
            self.db[dbfile]['t_max'] = max(teffs)

            # donor radius min/max
            self.db[dbfile]['r_min'] = min(radii)
            self.db[dbfile]['r_max'] = max(radii)

            self.db[dbfile]['total_time'] = ages[-1] - ages[0]

            self.logger.debug('trying to find the first onset of mass transfer for ' + dbfile)
            self.db[dbfile]['mt_onset'] = self.first_mt_start(mt_rates)

            self.logger.debug( 'm1_min, m1_max, m2_min, m2_max, mt_onset = '
                               + str(self.db[dbfile]['m1_min']) + ' ' + str(self.db[dbfile]['m1_max']) + ' '
                               + str(self.db[dbfile]['m2_min']) + ' ' + str(self.db[dbfile]['m2_max']) + ' '
                               + str(self.db[dbfile]['mt_onset'])  )

            if self.db[dbfile]['mt_onset'] < 0:
                nomt.append(dbfile)

        for dbfile in nomt:
            self.logger.info('pruning ' + dbfile + ' from the database for lack of MT')
            del self.db[dbfile]


    def get_vals( self, filepath: str ) -> dict:
        with open(filepath, 'rb') as infile:
            donor_masses, accretor_masses, mt_rates, periods, teffs, ages, radii, dts  = pickle.load(infile)
            return { 'm1': donor_masses
                     , 'm2': accretor_masses
                     , 'mt': mt_rates
                     , 'periods': periods
                     , 'teffs': teffs
                     , 'ages': ages
                     , 'radii': radii
                     , 'dts': dts }


    def view(self, query: ProgenitorQuery) -> dict:
        """
        Given a query, restrict the view of candidate data files to the ones within the mass and period range
        """

        query = query.query

        self.logger.info('creating a view for query: ' + str(query) )

        def viewfilter(dbfile: str, vals: dict, query: dict) -> bool:
            """Interval intersections for the requested query parameters
               and values recorded in the evolution track
            """

            self.logger.debug('comparing query parameters with db metadata for ' + dbfile )
            if vals['isbh'] != query['bhns']:
                self.logger.debug('compact object mismatch')
                return False
            if vals['m1_min'] > query['m1'][1]:
                self.logger.debug('min donor mass greater than max in the query')
                return False
            if vals['m1_max'] < query['m1'][0]:
                self.logger.debug('max donor mass smaller than the min in the query')
                return False
            if vals['m2_min'] > query['m2'][1]:
                self.logger.debug('min accretor mass grater than max in the query')
                return False
            if vals['m2_max'] < query['m2'][0]:
                self.logger.debug('max accretor mass smaller than min in the query')
                return False
            if vals['mt_min'] > query['mt'][1]:
                self.logger.debug('min MT rate greater than max in the query')
                return False
            if vals['mt_max'] < query['mt'][0]:
                self.logger.debug('max MT rate smaller than min in the query')
                return False
            if vals['p_min'] > query['p'][1]:
                self.logger.debug('min period greater than max in query')
                return False
            if vals['p_max'] < query['p'][0]:
                self.logger.debug('max period smaller than min in query')
                return False
            if vals['t_min'] > query['teff'][1]:
                self.logger.debug('min Teff exceeds max in query')
                return False
            if vals['t_max'] < query['teff'][0]:
                self.logger.debug('max Teff exceeds min in query')
                return False

            return True

        v = { dbfile:value
              for (dbfile, value)
              in self.db.items()
              if   viewfilter(dbfile, value, query) }

        self.logger.info('built view with ' + str(len(v)) + ' candidate data files' )
        self.logger.debug('the following files are included in the view:' + str(v.keys()))

        return v


    def first_mt_start(self, mt: list) -> int:
        """ Declare MT to have started if MT rate exceeds the threshhold over several
        consecutive values, defined by constants window and threshhold in the method"""

        window = 4
        threshhold = -15

        if len(mt) <= window + 1:
            self.logger.debug('track too small, assume MT starts at position 0')
            return 0

        for i in range( len(mt) - (window-1) ) :
            if min( mt[i:i+window] )  >= threshhold:
                self.logger.debug( 'MT starts at position ' + str(i) )
                return i

        self.logger.debug('no MT!')
        return -1


    def __str__(self) -> None:
        return str(self.db)


class ProgenitorSearch:

    progens = {}
    match_idx = {}

    def __init__( self
                 , query: ProgenitorQuery
                 , db: ProgenitorDatabase ) -> None:

        self.db = db
        self.query = query.query

        self.logger = logging.getLogger ('progentool')
        self.logger.info( 'starting search with query qry=' + str(self.query) )

        self.view =  db.view(query)
        self.do_search()


    def do_search(self):
        """Find files matching query parameters"""

        self.logger.debug('iterating over the '
                          + str( len(self.view.keys()) )
                          + ' candidate files in view')
        for dbfile in self.view.keys():
            vals = db.get_vals(dbfile)
            self.match_idx[dbfile] = [ idx for idx in range(len(vals['m1']))
                                       if match_props( vals['m1'][idx]
                                                       , vals['m2'][idx]
                                                       , vals['mt'][idx]
                                                       , vals['teff'][idx]
                                                       , vals['p'][idx]) ]

            self.progens[dbfile] = { 'm1_0': self.db[dbfile]['m1']
                                     , 'p_0': self.db[dbfile]['p']
                                     , 'm2_0': self.db[dbfile]['m2_min']
                                     , 'total_time': self.db[dbfile]['total_time']
                                     , 'observed_time': sum( [ 10**x
                                                               for x
                                                               in [ vals['dts'][idx]
                                                                    for idx in match_idx[dbfile] ] ] )
                                     , 'm1_at_mt_onset': vals['m1'][self.db[dbfile]['mt_onset']]
                                     , 'log10_p_at_mt_onset': np.log10( vals['p'][self.db[dbfile]['mt_onset']] ) }

    def match_props(data: dict) -> bool:
            query = self.query
            if (  data['m2']      < query['m2'][0]   or data['m2']   > query['m2'][1]
                  or data['mt']   < query['mt'][0]   or data['mt']   > query['mt'][1]
                  or data['p']    < query['p'][0]    or data['p']    > query['p'][1]
                  or data['teff'] < query['teff'][0] or data['teff'] > query['teff'][1]
                  or data['m1']   < query['m1'][0]   or data['m1']   > query['m1'][1]):
                return False
            else:
                return True

    def __str__(self):
        s = ''
        for progen in self.progens.values():
            for param in [ 'm1_0', 'p_0', 'm2_0', 'total_time', 'observed_time'
                           , 'm1_at_mt_onset', 'log10_p_at_mt_onset' ]:
                s += f'{progen[param]:.4f}' + ' '
            s += "\n"

        return (s)


if __name__ == "__main__":
    import sys
    import uuid

    req_id = str(uuid.uuid4())
    infile = sys.argv[1]
    db_location = os.environ['DB_LOCATION']

    logger = logging.getLogger('progentool')
    if 'DEBUG' in os.environ:
        logging.basicConfig( level = logging.DEBUG
                             , format = 'id=' +  req_id + ', t=%(asctime)s, message=%(message)s' )
    else:
        logging.basicConfig( level = logging.INFO
                             , format = 'id=' +  req_id + ', t=%(asctime)s, message=%(message)s' )

    logger.info('progenitor tool started. query input file: %s', infile )
    logger.info('database location: %s', db_location )

    db      = ProgenitorDatabase(db_location)
    query   = ProgenitorQuery(infile)
    results = ProgenitorSearch(query, db)

    print(results)
