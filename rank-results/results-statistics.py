import numpy as np
import glob

results_file = glob.glob("rank-results/complete list/results_*")

results_sts = []#np.empty([54,4])

color = ["blue", "red", "orange", "magenta", "green", "yellow", "brown", "black"]

for r in results_file:
	arr = np.loadtxt(r)
	r_algs = []
	for x in xrange(1, len(arr[1,:])):
		percent = (abs(arr[:,x]-arr[:,0]) == 0.);
		ratio = abs(arr[:,x]/arr[:,0]);
		r_algs.append([np.average(percent),np.std(percent),np.average(ratio),np.std(ratio)])
		pass
	
	results_sts.append(r_algs)
	pass

results_sts = np.array(results_sts)

print np.average(results_sts[:,1,0]), np.average(results_sts[:,2,0]), np.average(results_sts[:,3,0]), np.average(results_sts[:,4,0]), np.average(results_sts[:,5,0]), np.average(results_sts[:,6,0])

for x in xrange(0, len(results_sts[1,:])):
	string = "\plotbardata{" + color[x] + "}\n"
	string += "table [y] {\n"
	string += "x   y           error    label\n"
	string += "1   " + str(np.average(results_sts[:,x,0])) + " " + str(np.std(results_sts[:,x,0])) + " 2\n" 
	string += "};"

	print string
	pass

print "Ratio of approximation"
print "------------------------------------------"

for x in xrange(0, len(results_sts[1,:])):
	string = "\plotbardata{" + color[x] + "}\n"
	string += "table [y error=error] {\n"
	string += "x   y           error    label\n"
	string += "1   " + str(np.average(results_sts[:,x,2])) + " " + str(np.std(results_sts[:,x,2])) + " 2\n" 
	string += "};"

	print string
	pass