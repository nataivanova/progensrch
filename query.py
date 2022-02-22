import numpy as np
import pickle
import time
import os.path
import sys

# Progenitor Grid
m_arrl = np.arange(4.0,7.1,0.1)
m_arrs = np.arange(0.95,4.05,0.05)
p_arrl = np.arange(1.65,4.05,0.05)
p_arrs = np.arange(-0.60,1.66,0.02)

# Grid separation w.r.t density since simulations are stored in four folders
confs = {'lmlp':[m_arrl,p_arrl],'smlp':[m_arrs,p_arrl],'lmsp':[m_arrl,p_arrs],'smsp':[m_arrs,p_arrs]}

# m1 -> Donor Mass {Msol}
# m2 -> Accretor Mass {Msol}
# mt -> log10(MT Rate) {Msol/yr}
# p -> Orbital Period {days}
# teff -> log10(Donor Effective Temp) {K}
labels = ['m1','m2','mt','p','teff']

# Binary search algorithm used for Donor Mass
def search(array, element):
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

# Function to retrieve query information from .txt file and return a dictionary for the query
def get_query(qpath):
    
    try:
        qtemp = np.loadtxt(qpath,delimiter=",")
    except Exception as e:
        print("Error in query input, please recheck format")
        exit()

    #print(qtemp.shape)

    query = {}
    for i in range(qtemp.shape[0]):
        query[labels[i]] = qtemp[i]
    
    # Various input checks
    if query['m1'][0] < 0.0 or query['m1'][0] > query['m1'][1]:
        print("Wrong donor mass range")
        exit()
    if query['m2'][0] < 0.0 or query['m2'][0] > query['m2'][1]:
        print("Wrong BH mass range")
        exit()
    if query['mt'][0] < -50.0 or query['mt'][0] > query['mt'][1]:
        print("Wrong MT rate")
        exit()
    if query['p'][0] < 0.0 or query['p'][0] > query['p'][1]:
        print("Wrong orbital period")
    if query['teff'][0] > query['teff'][1]:
        print("Wrong Donor Teff")
        exit()

    return query

# Function to check whether a simulated system (data) matches all query properties
# Returns 1 if simulated system matches query and 0 if it doesn't
def match_props(data,query):
        
    if 'm1' in query:
        if data[0] < query['m1'][0] or data[0] > query['m1'][1]:
            #print("No m1")
            return 0
    if 'm2' in query:
        if data[1] < query['m2'][0] or data[1] > query['m2'][1]:
            #print("No m2")
            return 0
    if 'mt' in query:
        if data[2] < query['mt'][0] or data[2] > query['mt'][1]:
            #print("No mt")
            return 0
    if 'p' in query:
        if data[3] < query['p'][0] or data[3] > query['p'][1]:
            #print("No p")
            return 0
    if 'teff' in query:
        if data[4] < query['teff'][0] or data[4] > query['teff'][1]:
            #print("No teff")
            return 0
    return 1


# Function to look through data files to find progenitors for the query
def get_progens(query):
    
    progens = []
    
    # Loop through four folders
    for c in confs.keys():
        # Loops for initial mass and period values 
        for i in confs[c][0]:
            for j in confs[c][1]:
                
                # Ignore simulations where initial donor mass is less than the lower limit of donor mass in query
                if i >= query['m1'][0]:
                    try:
                        fpath = '/home/chatriks/pickle_runs_mesh/'+c+'/m_'+f'{i:4.2f}'+'_p_'+f'{j:4.2f}'+'.data'
                        if os.path.exists(fpath):
                            infile = open(fpath,'rb')
                            vals = pickle.load(infile)
                            infile.close()
                            
                            #print(fpath)
                            #print(vals.shape)
                            
                            # Ignore simulations where final donor mass is greater than upper limit of donor mass in query. Similarly, ignore simulations where initial accretor mass is greater than upper limit of accretor mass in query
                            if (vals[0][-1] <= query['m1'][1] and vals[1][0] <= query['m2'][1]):
                                flag = 0
                                # Indices in array to define window of entries where donor mass satisfies query donor mass limits
                                start_m = 0
                                end_m = 0
                                
                                # If upper limit of donor mass in query is greater than initial donor mass, start window at top of file
                                if vals[0][0] <= query['m1'][1]:
                                    start_m = 0
                                else:
                                # Binary search for upper donor mass limit
                                    start_m = search(vals[0],query['m1'][0])-1
                                
                                # If lower limit of donor mass in query is less than final donor mass end window at bottom of file
                                if vals[0][-1] >= query['m1'][0]:
                                    end_m = vals.shape[1]-1
                                else:
                                # Binary search for lower donor mass limit
                                    end_m = search(vals[0],query['m1'][0])
                                
                                # If donor mass window is found go through entries one by one to match all system properties to query. Flag all simulations with a matching system.
                                if (start_m <= end_m):
                                    #print("Mass matched ",start_m,end_m,end_m-start_m)
                                    for k in range(start_m,end_m+1):
                                        if (k < 0 or k >= vals.shape[1]):
                                            break
                                        flag = match_props([vals[0][k],vals[1][k],vals[2][k],vals[3][k],vals[4][k]],query)
                                        if flag == 1:
                                            progens.append([i,j])
                                            print("Found Progenitor -> "+fpath)
                                            break
                            
                        else:
                            print("No path found: "+fpath)

                    except Exception as e:
                        print("Error occurred in "+fpath+"\n"+str(e)+"\n")
    return progens

# Driver Code
print("Start Search")
start_time = time.time()

# Query .txt filename given as command line argument, change path according to preferences
qname = sys.argv[1]
qpath = ""+qname
query = get_query(qpath)
print("Query input taken:")
print(query)
progens = get_progens(query)
#print(progens)

# Progenitor list stored in corresponding .txt file, change path according to preferences
rpath = "progens_"+qname
outfile = open(rpath,'w')
for i in progens:
    s = f'{i[0]:4.2f}'+' '+f'{i[1]:4.2f}'+'\n' 
    outfile.write(s)
outfile.close()
print("Progenitor properties stored in: "+rpath)

print("Search Complete, time taken: "+str(time.time() - start_time)+" seconds")
