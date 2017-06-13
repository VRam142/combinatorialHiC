import os,sys,re
from collections import Counter

fhi = open(sys.argv[1])
fho = open(sys.argv[2], 'w') #File for easy plotting in ggplot

a = ['left', 'right', 'inner','outer']

long = {}
short = {}
dist = {"< 1000_human":0, "1000 - 20000_human":0, ">20000_human":0, "< 1000_mouse":0, "1000 - 20000_mouse":0, ">20000_mouse":0, "InterHuman":0, "InterMouse":0, "Interspecies":0}
for i in a:
	long[i] = 0
	short[i] = 0
total_reads = 0
for line in fhi:
	split = line.split()
	fcoord1 = int(split[1])
	total_reads +=1
	if fcoord1 == -1: continue
	rcoord2 = int(split[5])
	mapq1, mapq2 = int(split[7]), int(split[8])
	species1 = split[0]
	species2 = split[3]
	if species1 != species2:
		if re.search("human", species1):
			if re.search("human", species2):
				dist["InterHuman"] += 1
			elif re.search("mouse", species2):
				dist["Interspecies"] += 1
		else:
			dist["InterMouse"] += 1
		continue
	or1 = split[9]
	or2 = split[10]
	if rcoord2 - fcoord1 > 20000:
		if re.search("human", species1):
			dist[">20000_human"] += 1
		if re.search("mouse", species1):
			dist[">20000_mouse"] += 1
		if or1 == "-" and or2 == "-":
			long["left"] +=1
		if or1 == "+" and or2 == "+":
			long["right"] += 1
		if or1 == "+" and or2 == "-":
			long["inner"] += 1
		if or1 == "-" and or2 == "+":
			long["outer"] +=1
	else:
		if rcoord2 - fcoord1 < 1000:
			if re.search("human", species1):
				dist["< 1000_human"] += 1
			if re.search("mouse", species1):
				dist["< 1000_mouse"] += 1
                	if or1 == "-" and or2 == "-":
                        	short["left"] +=1
                	if or1 == "+" and or2 == "+":
                        	short["right"] +=1
                	if or1 == "+" and or2 == "-":
                        	short["inner"]+=1
                	if or1 == "-" and or2 == "+":
                        	short["outer"]+=1
		elif rcoord2 - fcoord1 >= 1000:
			if re.search("human", species1):
				dist["1000 - 20000_human"] +=1
			if re.search("mouse", species1):
				dist["1000 - 20000_mouse"] += 1
sum = 0
combined_species = Counter()
for i in dist:
	sum+= dist[i]
print "<H3>Summary Statistics for %s</H3>" % sys.argv[1]
print "<table border=\"4\">"
print "<tr><td>Filter</td><td>Distance_Class</td><td>Count</td><td>Percent</td><td>Total_Mapped</td><></tr>"
for i in dist:
	if i == "< 1000_human" or i == "< 1000_mouse":
		combined_species["< 1000"] += float(dist[i]) / float(sum)
	if i == ">20000_human" or i == ">20000_mouse":
		combined_species["> 20000"] += float(dist[i]) / float(sum)
	if i == "1000 - 20000_human" or i == "1000 - 20000_mouse":
		combined_species["1000 - 20000"] += float(dist[i]) / float(sum)
	if i == "InterHuman" or i == "InterMouse":
		combined_species["Interchromosomal"] += float(dist[i]) / float(sum)
	if i == "Interspecies":
		combined_species["Interspecies"] += float(dist[i]) / float(sum)
	print "<tr><td>%s</td><td>%s</td><td>%s</td><td>%.2f</td><td>%s</td></tr>" % ("Total",i,dist[i], 100 * float(dist[i]) / float(sum), sum)
print "</table>"
print ""
print "<p>Directionality breakdown</p>"
print "<table border=\"2\">"
print "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>" % ("Greater_than_20kb", long["left"], long["right"], long["inner"], long["outer"])
print "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>" % ("<= 1_kb", short["left"], short["right"], short["inner"], short["outer"])
print "</table>"
for i in combined_species:
	print >> fho, "%s\t%s" % (i, combined_species[i])

fhi.close()
fho.close()

