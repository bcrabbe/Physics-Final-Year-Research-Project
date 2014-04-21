#! /usr/bin/python
# -*- coding: utf-8 -*-
import fileinput as fi
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math as maths #localisation
import sys
import argparse

def normableHist(x,bins,normalise=False): #Because lxplus has a really old numpy so hist lacks the density option
	h, b = np.histogram(x,bins);
	if normalise:
		s = float(h.sum())
		#print s
		h = [num/s for num in h]
	return h, b

def cleanList(l):
	r = []
	for item in l:
		try:
			#floats = [float(string) for string in item.split()]
			#f = floats[3]
			#f = floats[0]/floats[1]
			f  = float(item)
			if (not maths.isnan(f)) and (args.ll < f  < args.ul):
				r.append(f)
		except:
			continue
	return r

parser = argparse.ArgumentParser(description='Plot data in hitograms');

parser.add_argument('-b',help='number of bins (defaults to 500)',type=int,default= 500)

parser.add_argument('-n',help='normalise frequencies and plot all on the same figure (defaults to False)',type=bool,nargs='?',default=False,const=True);

parser.add_argument('-l',help='Make the y scale(s) logarithmic (defaults to False)',type=bool,nargs='?',default=False,const=True)
parser.add_argument('-c',help='add a cut',type = float)

graphTypeGroup = parser.add_mutually_exclusive_group(required=True);
graphTypeGroup.add_argument('-B',help='draw as a bar graph',action='store_true',default=False)
graphTypeGroup.add_argument('-L',help='draw as a line graph',action='store_true',default=False)

parser.add_argument('-i',help='input file(s) (defaults to stdin)',type=argparse.FileType('r'),default=[sys.stdin],nargs='*')
parser.add_argument('-ll',help='lower limit',type=float,default=-sys.float_info.max)
parser.add_argument('-ul',help='upper limit',type=float,default=sys.float_info.max)
args = parser.parse_args()
print args
matplotlib.rcParams['text.latex.unicode']=True #for greek letters

for inFile in args.i:
	cutPassedProportions = []
	arrayList = cleanList( inFile.readlines()) #return a list of ~nan floats
	numItems = len(arrayList)
	numFailedItems = 0
	array = np.array(arrayList)
	if args.c != None:
		for item in array:
			if item < float(args.c):
				numFailedItems += 1
				labelStr = inFile.name.decode('utf-8') + ': ' + str(float(100*numFailedItems)/float(numItems)) + '%'
	else:
		labelStr = inFile.name.decode('utf-8')
	#print inFile.name.decode('utf-8') + ': ' + str(numFailedItems) + ','+ str(numItems)
	hist,bins = normableHist(array,bins=args.b,normalise=args.n)
	width = bins[1] - bins[0]
	if not args.n:
		plt.figure()
		plt.title(inFile.name)
	#print hist.length, bins.length
	s = 0
	for num in hist:
		s += num
	if args.l:
		plt.yscale('log')
	if args.L:
		plt.plot(bins[:-1], hist,label=labelStr)
	if args.B:
		plt.bar(bins[:-1],height=hist,label=labelStr,width=bins[1]-bins[0])
	plt.xlabel("E/P ratio");
	plt.ylabel("Relative frequency");
	
		
if args.c != None:
	plt.axvline(args.c,label='cut',color='black')
	
plt.legend(loc=0)
plt.show()
