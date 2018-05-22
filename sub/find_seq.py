# A script to shuffle the sequence of genes to find the optimal plot

#############################################################################
# A tool for finding the optimal seq. of the genes on xaxis                 #       
#############################################################################
# authors:                                                                  #
# Elmer stickel                                                             #
# Netherlands Cancer Insitute                                               #
# 8th March 2016                                                            #
# Contact: elmer.stickel@posteo.net                                         #
# Commercial use of the software is prohibited without                      #
# prior approval by Elmer stickel                                           #
#############################################################################
# The Python implementation of the Knuth Durstenfeld algoritm  used in this #
# script was obtained from the GitHub repository of Joe Di Castro.          #
# See https://github.com/joedicastro/python-recipes/blob/master/shuffle.py  #
#############################################################################

import math
import random
import os
import pandas as pd
import numpy as np
import sys
import timeit

start = timeit.default_timer()

n=25
graph_width=1200
graph_heigh=800
maxy=375 # This is the max p-value (in negative number to where the y-axis should reach, take a bit extra))
limit=25000 # The upper limit for the number of iterations to prevent the program to run for infite amount of time if overlap cannot be resolved 

# Read the CSV-file
infile=sys.argv[1]
outfile=sys.argv[2]
csvfile=pd.read_csv(infile, sep="\t", header=None)
csvfile = csvfile.reset_index()
maxins=csvfile[3].max() # Find max value in number of insertions required for scaling of the bubbles
ipd=(0.05*graph_width)/(np.sqrt(maxins/math.pi))
csvfile['radius']=np.sqrt(csvfile[3]/math.pi)*ipd
# Get length of dataframe
length = len(csvfile.index)
# Extract the id and p-value from the csv-file and put them in a Pandas Series
#rawdata=DataFrame(data=(csvfile.iloc[:,[3,6]]), index=idx)
rawdata=csvfile.iloc[:,[0,7,8]]
rawdata.columns=['index', 'fcpv', 'radius']
# Sort the Pandas Series on p-values
sorteddata=rawdata.sort_values(by='fcpv')
# Extract the top n-hits (as defined by 'n'), yield a new Pandas Series
hits=sorteddata[:n]
dfhits=pd.DataFrame(hits)
# fixlst holds the indexes of the sorted list, so when matching these to their new indexes we can find back the corresponding gene
fixlst=sorteddata.index.tolist()
# lst holds the new indexes assigned to the genes by the shuffle function
lst=sorteddata.index.tolist()
# Pixel density x-axis and y-axis
pdx=graph_width/float(length)
pdy=graph_heigh/maxy

def knuth_durstenfeld(lst):
    for idx in reversed(range(1, len(lst))):
        # pick an element in lst[:idx+1] with which to exchange lst[idx]
        sel = random.randrange(idx + 1)
        lst[idx], lst[sel] = lst[sel], lst[idx]
        

# Execute the shuffle function on lst
best=-200
count=0
while (best<0 and count<limit):
    knuth_durstenfeld(lst)
    tmp=pd.concat([pd.Series(fixlst), pd.Series(lst)], axis=1)
    tmp.columns=['index', 'seq']
    dftmphits=tmp[:n]
    tmp.columns=['index', 'seq']
    cc=pd.merge(dftmphits, dfhits, on='index')
    cc['x']=cc['seq']*float(pdx)
    cc['y']=-np.log10(cc['fcpv'])*pdy
    scalars = []
    a=0
    b=0
    s=pd.Series()
    for a in range(0, n):
        ar=cc.iat[a,3]
        ax=cc.iat[a,4]
        ay=cc.iat[a,5]
        for b in range(a+1, n):
            br=cc.iat[b,3]
            bx=cc.iat[b,4]
            by=cc.iat[b,5]
            # There's no need to make absolute values because when we put them in Pythagoras, the power of 2 cancels out the effect of the minus
            dx=ax-bx # delta x
            dy=ay-by # delta y
            dxy=math.sqrt(dx*dx+dy*dy)
            cdxy=dxy-ar-br
            scalars.append(cdxy)
    s=pd.Series(scalars)
    if (s.min()>best):
        best=s.min()
        bestseq=tmp.copy()
    count += 1
finaldf=pd.merge(csvfile, bestseq, on='index')
finaldf=finaldf.drop('index', 1)
finaldf.to_csv(outfile, sep="\t", header=None, index=False)
stop = timeit.default_timer()
speed = int(round(count/(stop - start)))
if (count < limit):
    print("Found optimal seqeuence in %i iterations @ %i iterations/sec." %(count,speed))
else:
    print("Reached max number of %i iterations. The nearest optimal seqeuence yields %spx overlap. Speed: %i iterations/sec." %(limit,best,speed))
