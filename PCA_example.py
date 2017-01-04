#Filtered human cells, PCA on the raw interchromosomal matrices (logged counts, ignoring the diagonal)
%matplotlib inline
import numpy
import matplotlib.pyplot as plt
from scipy.stats import spearmanr


ML1_coverage = open(sys.argv[1]) #Valid coverage file ML1
ML2_coverage = open(sys.argv[2]) #Valid coverage file ML2
PL1_coverage = open(sys.argv[3]) #Valid coverage file PL1
PL2_coverage = open(sys.argv[4]) #Valid coverage file PL2

coverages = [ML1_coverage,ML2_coverage,PL1_coverage,PL2_coverage]

labels = []
covs = []
replicate = []
repid = 0
for i in coverages:
    repid += 1
    for line in i:
        replicate.append(repid)
        cov = int(line.split()[0])
        species = line.split()[1]
        if species == "HeLa":
            labels.append("blue")
        if species == "HAP1":
            labels.append("yellow")
        covs.append(cov)
    i.close()

ML1_matrix = numpy.loadtxt(sys.argv[5], delimiter=",") #ML1 matrix
ML2_matrix = numpy.loadtxt(sys.argv[6], delimiter=",") #ML2 matrix
PL1_matrix = numpy.loadtxt(sys.argv[7], delimiter=",") #PL1 matrix
PL2_matrix = numpy.loadtxt(sys.argv[8], delimiter=",") #PL2 matrix

total = np.vstack((ML1_matrix, ML2_matrix, PL1_matrix, PL2_matrix))

pca = PCA(n_components=8)
PCA_test = pca.fit_transform(total)

plt.scatter(PCA_test[:,0], covs, s = 60, c=labels)
plt.show()

plt.scatter(PCA_test[:,0], PCA_test[:,1], s = 10, c=labels)
plt.show()

plt.scatter(PCA_test[:,2], PCA_test[:,3], s = 10, c=labels)
plt.show()

for i in range(len(PCA_test[:,0])):
    print "%s\t%s\t%s\t%s\t%s\t%s" % (PCA_test[:,0][i], PCA_test[:,1][i], PCA_test[:,2][i], labels[i], covs[i], replicate[i])




