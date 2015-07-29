import numpy as np
import glob

results_file = glob.glob("rank-results/time/results_*")

results_sts = []#np.empty([54,4])

color = ["blue", "red", "orange", "magenta", "green", "yellow", "brown", "black"]

for r in results_file:
	arr = np.loadtxt(r)
	r_algs = []

	saida = str()
	saida = str(np.average(arr[:,0])) 
	for x in [1,2,5,7,9]:
		saida += " & " + "{:.1e}".format(np.average(arr[:,x])) + "$\\pm$(" + "{:.1e}".format(np.std(arr[:,x])) + ")"
		pass
	#print saida
	results_sts.append(r_algs)
	pass

for r in results_file:
	arr = np.loadtxt(r)
	r_algs = []

	saida = str()
	saida = str(np.average(arr[:,0])) 
	for x in xrange(3,10):
		saida += "  " + "{:.1e}".format(np.average(arr[:,x])) + " " + "{:.1e}".format(np.std(arr[:,x]))
		pass
	print saida
	results_sts.append(r_algs)
	pass