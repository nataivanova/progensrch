#!/usr/bin/env python3

import numpy as np
import os
import pickle

class ProgenSearchException(Exception):
    errors = []

    def __init__(self, errors):
        self.errors = errors

    def __str__(self):
        return str(self.errors)

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
    labels = ['m1','m2','mt','p','teff']

    query = {}

    def __init__(self, qpath):
        print (qpath)
        try:
            q = open(qpath,"r")
            bhns = (q.readline()).split(None, 1)[0]
            q.close()

            qtemp = np.loadtxt(qpath,delimiter=",",skiprows=1)
        except Exception as e:
            print("Error in query input, please recheck format")
            print(str(e))
            exit()

        #print(qtemp.shape)
        self.query['bhns'] = int(bhns)
        for i in range(qtemp.shape[0]):
            self.query[self.labels[i]] = qtemp[i]
        
        return None


class ProgenitorSearcher:

    query = {}
    confs = {}
    db_location = '/tmp/'
    progens = []
    errors  = []

    def __init__(self, query, db_location):

        self.db_location = db_location

        self.query = query
        self.validate(self.query)

        # Progenitor Grid
        m_arrl = np.arange(4.0,7.1,0.1)
        m_arrs = np.arange(0.95,4.05,0.05)
        p_arrl = np.arange(1.65,4.05,0.05)
        p_arrs = np.arange(-0.60,1.66,0.02)

        self.bh_masses = [7.0,10.0]

        # Grid separation w.r.t density since simulations are stored in four folders
        self.confs = { 'lmlp':[m_arrl,p_arrl]
                     , 'smlp':[m_arrs,p_arrl]
                     , 'lmsp':[m_arrl,p_arrs]
                     , 'smsp':[m_arrs,p_arrs] }


    def validate(self, query):
        # Various input checks
        if (query['m1'][0] < 0.0 or query['m1'][0] > query['m1'][1]):
            self.errors.append("Wrong donor mass range")
        if (query['m2'][0] < 0.0 or query['m2'][0] > query['m2'][1]):
            self.errors.append("Wrong BH mass range")
        if (query['mt'][0] < -50.0 or query['mt'][0] > query['mt'][1]):
            self.errors.append("Wrong MT rate")
        if (query['p'][0] < 0.0 or query['p'][0] > query['p'][1]):
            self.errors.append("Wrong orbital period")
        if query['teff'][0] > query['teff'][1]:
            self.errors.append("Wrong Donor Teff")

        if (self.errors):
            raise QueryParametersException(self.errors)

    def match_props(self, data):
        query = self.query
        # print (query)
        try:
            if 'm1' in query:
                if (data[0] < query['m1'][0] or data[0] > query['m1'][1]):
                    print("No m1")
                    return 0
            if 'm2' in query:
                if (data[1] < query['m2'][0] or data[1] > query['m2'][1]):
                    print("No m2")
                    return 0
            if 'mt' in query:
                if (data[2] < query['mt'][0] or data[2] > query['mt'][1]):
                    print("No mt")
                    return 0
            if 'p' in query:
                if (data[3] < query['p'][0] or data[3] > query['p'][1]):
                    print("No p")
                    return 0
            if 'teff' in query:
                if (data[4] < query['teff'][0] or data[4] > query['teff'][1]):
                    print("No teff")
                    return 0
        except Exception as e:
            print("Error in matching properties")
            return 0

        return 1

    # Binary search algorithm used for Donor Mass
    def search(self, array, element):
        
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


    # Function to find index where system starts MT
    def find_mt_start(find, mt_arr):
        
        idx_start = 0

        if len(mt_arr) <= 5:
            return idx_start

        for i in range(len(mt_arr)-3):
            if (mt_arr[i] >= -15 and mt_arr[i+1] >= -15 and mt_arr[i+2] >= -15 and mt_arr[i+3] >= -15):
                idx_start = i
                break
        
        return idx_start

    def gen_search_paths(self):
        [ print (  '/srv/progen_tool/'
                 + str(self.query['bhns'])    + '/'
                 + str(accretor_mass)         + '/'
                 + str(quadrant)         + '/'
                 + str(period)           + '/'
                 + str(donor_mass) )
          for accretor_mass in [ 5, 7, 10 ]
          for quadrant   in self.confs.keys()
          for donor_mass in self.confs[quadrant][0] if donor_mass >= self.query['m1'][0]
          for period     in self.confs[quadrant][1]
          ]

    # Function to look through data files to find progenitors for the query
    def do_search(self):
        query = self.query
        confs = self.confs
        
        progens = []
        
        data_paths = []
        # Select relevant databases to search through
        if 0 == query['bhns']:
            data_paths.append(self.db_location)
        else:
            for x in self.bh_masses:
                if(x > query['m2'][1]):
                    continue
                if(x < query['m2'][0]-3.5):
                    continue
                data_paths.append(self.db_location + '/runs' +str(int(x))+'_data/')
        print(data_paths)

        print(confs)
        [print(x) for x in confs.keys()]

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
                                    if (vals[0][-1] <= query['m1'][1] and vals[1][0] <= query['m2'][1]):
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
                                            start_m = self.search(vals[0],query['m1'][0])-1
                                        
                                        # If lower limit of donor mass in query is less than final donor mass end window at bottom of file
                                        if vals[0][-1] >= query['m1'][0]:
                                            end_m = vals.shape[1]-1
                                        else:
                                        # Binary search for lower donor mass limit
                                            end_m = self.search(vals[0],query['m1'][0])
                                        
                                        # If donor mass window is found, go through simulation entries one by one to match all system properties to query. Flag all simulations with a matching system
                                        if (start_m <= end_m):
                                            #print("Mass matched ",start_m,end_m,end_m-start_m)
                                            idx_low = -1
                                            idx_high = -1
                                            obs_time = 0.0
                                            mt_start = 0

                                            for k in range(start_m,end_m+1):
                                                if (k < 0 or k >= vals.shape[1]):
                                                    break
                                                flag = self.match_props([vals[0][k],vals[1][k],vals[2][k],vals[3][k],vals[4][k]])
                                                
                                                # If a match is found, find the time spent as observed system, each traversal across the window is added to total observed time incrementally
                                                if (flag == 1 and idx_low == -1):
                                                    idx_low = k

                                                if ((flag == 0 and idx_low != -1 and idx_high == -1) or (flag == 1 and k == end_m and idx_high == -1)):
                                                    pflag = 1
                                                    idx_high = k
                                                    obs_time = obs_time + vals[5][idx_high]-vals[5][idx_low]
                                                    idx_low = -1
                                                    idx_high = -1

                                            # If simulated system spends time as the observed system, find start of MT information and add to list of progenitors
                                            if pflag == 1:
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
                    + f'{i[4]:.4f}' +' '
                    + f'{i[5]:.4f}' +' '
                    + f'{i[6]:.4f}'
                    +'\n')

        return (s)

if __name__ == "__main__":
    import sys
    import uuid
    x = ProgenitorQuery(sys.argv[1])
    db_location = os.environ['DB_LOCATION']
    try:
        s = ProgenitorSearcher(x.query, db_location)
    except QueryParametersException as e:
        print(e)
        exit (1)

    s.gen_search_paths()
    try:
        s.do_search()
    except Exception as e:
        print (type(e), str(e))
    print(s)
    req_id = str(uuid.uuid4())
    outpath = '/tmp/' + req_id + '.result'
    outfile = open(outpath, 'wb')
    outfile.write(bytes(str(s), 'UTF-8'))
    outfile.flush()
    outfile.close()
    outfile = open(outpath, 'rb')
    for i in outfile:
        print(i)
