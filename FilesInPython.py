#!/usr/bin/python

import sys
import time
import csv
import math
import random
from bisect import insort
import scipy
from scipy.stats import histogram
from matplotlib import pyplot


#_______________________________________________________________________________
def write_data(file_name):
    '''
        Writes some data to a csv data file

        *** You may want to look at the csv.DictWriter module which does
            basically the same this as this function.
    '''

    # Generate some fake data
    print "generating data..."
    data = []
    keys = ['type_a','type_b','type_c']
    for row in range(10000):
        idx = random.randint(0,len(keys)-1)
        val = random.gauss(0,1)
        data.append((keys[idx], val))
    
    # Write the fake data to a csv file
    print "opening file..."
    f = open(file_name,"wb")
    print "writing data..."

    # Start timer
    t1 = time.time()

    # Load data (but keep it sorted)
    writer = csv.writer(f, delimiter=',')
    for n in data:
        writer.writerow(n)
    
    # Stop timers
    t2 = time.time()
    f.close()

    print "time taken:", int((t2 - t1)*10)/10.0, "sec"

#_______________________________________________________________________________
def read_data(file_name, do_sort=False):
    '''
        Takes a csv reader and reads in the data into a dictionary,
        where the keys.

        *** You may want to look at the csv.DictReader module which does
            basically the same this as this function.
    '''

    # Read data from a csv file
    print "opening file..."
    f = open(file_name,"rb")
    print "loading data..."

    # Start timer
    t1 = time.time()

    # Load data (but keep it sorted)
    reader = csv.reader(f, delimiter=',')
    dataDict = {}
    for row in reader:

        # Data file format dependent code
        key = row[0]
        val  = float(row[1])

        if do_sort:
            insort(dataDict.setdefault(key, []), val)
        else:
            dataDict.setdefault(key, []).append(val)

    # Stop timers
    t2 = time.time()
    f.close()

    print "loaded", reader.line_num, "rows,"
    print "time taken:", int((t2 - t1)*10)/10.0, "sec"

    return dataDict


#_______________________________________________________________________________
def rm_outliers(data, lq, uq, is_sorted=False):
    '''
        Removes the outliers from a set of data.

        Inputs: data    The data list
                lq      The lower quartile percentage (should be less than 1)
                uq      The upper quartile percentage (should be less than 1)
    '''
    if not is_sorted: data.sort()

    n = float(len(data))
    lbidx = int(round(n*lq))
    ubidx = int(round(n*uq))
    
    if lbidx >= 0:
        lb = data[lbidx]
    else:
        lb = data[0]
        
    if ubidx < len(data):
        ub = data[ubidx]
    else:
        ub = data[len(data) - 1]

    while data[0] < lb:
        data.pop(0)

    i = len(data) - 1
    while data[i] > ub:
        data.pop(i)
        i -= 1

#_______________________________________________________________________________
def get_stats(data, is_sorted=False):
    '''
        Returns a dictionary of statistics for a list of numerical data
    '''
    stats = {}
    if not is_sorted: data.sort()
    n = len(data)
    stats['count'] = n
    stats['min']   = data[0]
    stats['max']   = data[n - 1]
    stats['mean'] = scipy.mean(data)
    stats['stdv'] = scipy.std(data, ddof=1)

    return stats

#_______________________________________________________________________________
# Main function
def main(argv):
     # Get input file name
    try:
        file_name = sys.argv[1]

    except:
        print "USAGE: python FilesInPython.py [file name]"
        sys.exit(-1)

    # Make a fake data file
    write_data(file_name)

    # Read data into a dictionary
    dataDict = read_data(file_name)

    # Remove outliers
    print "removing outliers..."
    for key in dataDict:
        rm_outliers(dataDict[key], 0.01, 0.99)

    # Calculate statistics for each key
    for key in dataDict:
        print "--------------------------------------"
        print key+":"

        # Print statistics
        stats = get_stats(dataDict[key])
        statList = sorted(stats.keys())
        for stat in statList:
            print stat+":  "+str(stats[stat])

        # Display histogram
        nbins = 10
        fig = pyplot.figure(figsize=(8, 6),
                            dpi=100,
                            facecolor='w',
                            edgecolor='k')
        pyplot.hist(dataDict[key], bins=nbins, color='b')
        pyplot.ylabel("count",figure=fig)
        pyplot.title(key,figure=fig)
        pyplot.show()


#_______________________________________________________________________________
if __name__ == '__main__':
    main(sys.argv)

