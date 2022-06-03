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
            self.init_cache(cache_location)
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

    def init_cache(self):
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

            self.logger.debug('trying to find the first onset of mass transfer for ' + dbfile)
            self.db[dbfile]['mt_onset'] = self.first_mt_start(mt_rates)

            self.logger.debug( 'm1_min, m1_max, m2_min, m2_max, mt_onset = '
                               + str(self.db[dbfile]['m1_min']) + ' ' + str(self.db[dbfile]['m1_max']) + ' '
                               + str(self.db[dbfile]['m2_min']) + ' ' + str(self.db[dbfile]['m2_max']) + ' '
                               + str(self.db[dbfile]['mt_onset'])  )

            if self.db[dbfile]['mt_onset'] < 0:
                nomt.append(dbfile)

        for dbfile in nomt:
            self.logger.info('removing ' + dbfile + ' from the database for lack of MT')
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

    progens = []

    def __init__(self, query: ProgenitorQuery, db: ProgenitorDatabase) -> None:

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
                          +  ' candidate files in view')
        for dbfile in self.view.keys():
            pass


    # Binary search algorithm in a list pre-sorted in ascending order
    def search(self, array, element) -> int:

        if element <= array[0]:
            return 0
        if element >= array[-1]:
            return len(array) - 1

        mid = 0
        start = 0
        end = len(array)
        step = 0

        while (start <= end):
            step = step+1
            mid = (start + end) // 2

            if mid == len(array):
                return len(array)-1

            if element == array[mid]:
                return mid

            if element > array[mid]:
                end = mid - 1
            else:
                start = mid + 1
        return start



    def extract_progens(self, vals) -> None:

        flag, pflag, start_m1, end_m1 = 0, 0, 0, 0
        start_m1, end_m1 = self.search(vals['m1'], query['m1'][0]),   self.search(vals['m1'], query['m1'][1])

        idx_low = -1
        idx_high = -1
        obs_time = 0.0
        mt_start = 0
        for k in range(start_m, end_m+1):
            if match_props( [ vals[_][k] for _ in ['m1', 'm2', 'mt', 'p', 'teff'] ] ):
                idx_low = k
            if (   (flag    and idx_low != -1 and idx_high == -1)
                   or (flag and k == end_m    and idx_high == -1) ):
                pflag = 1
                idx_high = k
                obs_time = obs_time + vals[5][idx_high]-vals[5][idx_low]
                idx_low = -1
                idx_high = -1
                # If simulated system spends time as the observed system, find start of MT information and add to list of progenitors
            if pflag:
                mt_start = self.find_mt_start(vals[2])
                tot_time = vals[5][-1]-vals[5][0]
                progens.append([i,j,vals[1][0],obs_time,tot_time,vals[0][mt_start],np.log10(vals[3][mt_start])])
                print("Found Progenitor -> "+fpath)
            else:
                print("No path found: "+fpath)

        self.progens = progens
        if not progens:
            raise NoProgensFoundException('No progenitors found')


    def delme(self):
        progens = []

        # Loop through databases
        for dpath in data_paths:
            # Loop through four folders
            for c in confs.keys():
                # Loops for initial mass and period values
                for i in confs[c][0]:
                    for j in confs[c][1]:

                        # Ignore simulations where initial donor mass is less than the lower limit of donor mass in query
                        if i >= query['m1'][0]:
                            try:
                                fpath = dpath+c+'/m_'+f'{i:4.2f}'+'_p_'+f'{j:4.2f}'+'.data'
                                if os.path.exists(fpath):
                                    infile = open(fpath,'rb')
                                    vals = pickle.load(infile)
                                    infile.close()

                                    #print(fpath)
                                    #print(vals.shape)

                                    # Ignore simulations where final donor mass is greater than upper limit of donor mass in query. Similarly, ignore simulations where initial accretor mass is greater than upper limit of accretor mass in query
                                    if (vals[0][-1] <= query['m1'][1] and vals[1][-1] >= query['m2'][0]):
                                        # Flags to indicate if system is a progenitor
                                        flag = 0
                                        pflag = 0
                                        # Indices in array to define window of entries where donor mass satisfies query donor mass limits
                                        start_m = 0
                                        end_m = 0

                                        # If upper limit of donor mass in query is greater than initial donor mass, start window at top of file
                                        if vals[0][0] <= query['m1'][1]:
                                            start_m = 0
                                        else:
                                        # Binary search for upper donor mass limit
                                            start_m = self.search(vals[0],query['m1'][1])-1

                                        # If lower limit of donor mass in query is less than final donor mass end window at bottom of file
                                        if vals[0][-1] >= query['m1'][0]:
                                            end_m = vals.shape[1]-1
                                        else:
                                        # Binary search for lower donor mass limit
                                            end_m = self.search(vals[0],query['m1'][0])

                                        # If donor mass window is found, go through simulation entries one by one to match all system properties to query. Flag all simulations with a matching system
                                        if (start_m <= end_m):
                                            #print("Mass matched ",start_m,end_m,end_m-start_m)
                                            obs_time = 0.0
                                            mt_start = 0

                                            for k in range(start_m,end_m+1):
                                                if (k < 0 or k >= vals.shape[1]):
                                                    break
                                                flag = match_props([vals[0][k],vals[1][k],vals[2][k],vals[3][k],vals[4][k]])

                                                # If a match is found, find the time spent as observed system, each traversal across the window is added to total observed time incrementally
                                                if flag:
                                                    pflag = 1
                                                    obs_time += 10**vals[7][k]

                                            # If simulated system spends time as the observed system, find start of MT information and add to list of progenitors
                                            if pflag:
                                                mt_start = self.find_mt_start(vals[2])
                                                tot_time = vals[5][-1]-vals[5][0]

                                                progens.append([i,j,vals[1][0],obs_time,tot_time,vals[0][mt_start],np.log10(vals[3][mt_start])])
                                                print("Found Progenitor -> "+fpath)

                                else:
                                    print("No path found: "+fpath)

                            except Exception as e:
                                print("Error occurred in "+fpath+"\n"+str(e)+"\n")

        self.progens = progens
        if not progens:
            raise NoProgensFoundException('No progenitors found')


    def __str__(self):
        s = ''
        for i in self.progens:
            s += (    f'{i[0]:.2f}' + ' '
                    + f'{i[1]:.2f}' + ' '
                    + f'{i[2]:.2f}' + ' '
                    + f'{i[3]:.4f}' + ' '
                    + f'{i[4]:.4f}' + ' '
                    + f'{i[5]:.4f}' + ' '
                    + f'{i[6]:.4f}'
                    +'\n')

        return (s)


if __name__ == "__main__":
    import sys
    import uuid

    req_id = str(uuid.uuid4())
    infile = sys.argv[1]
    db_location = os.environ['DB_LOCATION']

    logger = logging.getLogger('progentool')
    logging.basicConfig( level = logging.INFO
                        , format = 'id=' +  req_id + ', t=%(asctime)s, message=%(message)s' )
    logger.info('progenitor tool started. query input file: %s', infile )
    logger.info('database location: %s', db_location )

    db      = ProgenitorDatabase(db_location)
    query   = ProgenitorQuery(infile)
    results = ProgenitorSearch(query, db)

    exit(0)

    try:
        s = ProgenitorSearcher(qry.query, db_location)
    except QueryParametersException as e:
        logger.error(e)
        exit (1)

    try:
        s.do_search()
    except Exception as e:
        print (type(e), str(e))
    print(s)
    outpath = '/tmp/' + req_id + '.result'
    outfile = open(outpath, 'wb')
    outfile.write(bytes(str(s), 'UTF-8'))
    outfile.flush()
    outfile.close()
    outfile = open(outpath, 'rb')
    for i in outfile:
        print(i)
