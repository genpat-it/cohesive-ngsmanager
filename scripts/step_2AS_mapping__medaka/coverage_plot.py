#!/usr/bin/env python3
import sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl 

mpl.rcParams['agg.path.chunksize'] = 10000

def plotCov(depthFile,plotFile):
	x = []
	y = []

	table = open(depthFile, 'r')
	for tline in table:
		line = tline.split("\t")
		x.append(int(line[1]))
		y.append(int(line[2]))
	avg = np.average(y)
	plt.plot(x, y, 'k')
	plt.xlabel('Position in Genome')
	plt.ylabel('Depth of Coverage')
	plt.hlines(avg, 0, 30000, color="red", linestyles='solid')
	plt.text(30000,avg, "AVG depth",color="red",ha="left", va="center")
	# plt.show()
	plt.savefig(plotFile,dpi= 1000)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.stderr.write(
            "Usage: %s <samtools_depth_out> <plot_name>" % (sys.argv[0]))
        sys.exit(1)
    else:
        plotCov(sys.argv[1], sys.argv[2])