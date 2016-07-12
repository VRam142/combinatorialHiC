import os,sys,re
from collections import Counter
import random

fhi = open(sys.argv[1]) #Open filehandle of associations bedpe
cnt = Counter() #set up a Counter object for counting the occurence of barcode pairs
frag_occurrence = {}
bc = {} #Dict for barcode associatins
bc_long = {} #that are above a certain distance class
bc_random = {} #and for randoms as
bc_random_long = {} #as well
bcs = [] #Store all the barcode labels so you can shuffle them later

for line in fhi: #Create a list of barcode combinations and their associated counts for later
	split = line.split()
	cnt["%s-%s"  % (split[11], split[12])] += 1
	bcs.append("%s-%s" % (split[11], split[12]))

random.shuffle(bcs)



for combo in cnt:
	frag_occurrence[combo] = Counter()
	bc[combo] = {}
	bc_random[combo] = {}
	bc_long[combo] = {}
	bc_random_long[combo] = {}
	bc[combo]["human"] = 0
	bc[combo]["mouse"] = 0
	bc_long[combo]["human"] = 0
	bc_long[combo]["mouse"] = 0
	bc_random[combo]["human"] = 0
	bc_random[combo]["mouse"] = 0
	bc_random_long[combo]["human"] = 0
	bc_random_long[combo]["mouse"] = 0
fhi.close()

fhi = open(sys.argv[1]) #Reopen the filehandle
i = 0 #Counter for random barcodes
for line in fhi:
	split = line.split()
	barcode = "%s-%s" % (split[11], split[12])
	fragid1 = split[13]
	fragid2 = split[15]
	strand1, strand2 = split[9:11]
	dist1, dist2 = split[14], split[16]
	if fragid1 == fragid2:
		frag_occurrence[barcode][(fragid1, strand1, dist1)] += 1
	elif split[1] == split[4]:
		frag_occurrence[barcode][(fragid1, strand1, dist1)] += 1	
	else:
		frag_occurrence[barcode][(fragid1, strand1, dist1)] += 1
		frag_occurrence[barcode][(fragid2, strand2, dist2)] += 1
	random_bc = bcs[i]
	if re.search("human", split[0]):
		if re.search("human", split[3]):
			if split[0] == split[3]:
				bc[barcode]["human"] += 1
				bc_random[random_bc]["human"] += 1
				if int(split[4]) - int(split[1]) > 1000:
					bc_long[barcode]["human"] += 1
					bc_random_long[random_bc]["human"] += 1
			elif split[0] != split[3]:
				bc[barcode]["human"] += 1
				bc_random[random_bc]["human"] += 1
				bc_long[barcode]["human"] += 1
				bc_random_long[random_bc]["human"] += 1
	elif re.search("mouse", split[0]):
		if re.search("mouse", split[3]):
			if split[0] == split[3]:
				bc[barcode]["mouse"] += 1
				bc_random[random_bc]["mouse"] += 1
				if int(split[4]) - int(split[1]) > 1000:
					bc_long[barcode]["mouse"] += 1
					bc_random_long[random_bc]["mouse"] += 1
			elif split[0] != split[3]:
				bc[barcode]["mouse"] += 1
				bc_random[random_bc]["mouse"] += 1
				bc_long[barcode]["mouse"] += 1
				bc_random_long[random_bc]["mouse"] += 1
	i += 1
total = 0
for combo in cnt:
	total = bc[combo]["human"] + bc[combo]["mouse"]
	total_random = bc_random[combo]["human"] + bc_random[combo]["mouse"]
	if total == 0: continue
	frac_human = float(bc[combo]["human"]) / float(total)
	frac_mouse = float(bc[combo]["mouse"]) / float(total)
	if total_random == 0:
		frac_human_random = 0
		frac_mouse_random = 0
	else:
		frac_human_random = float(bc_random[combo]["human"]) / float(total_random)
		frac_mouse_random = float(bc_random[combo]["mouse"]) / float(total_random)
	bc1, bc2 = combo.split("-")
	ones = 0
	twos = 0
	threes = 0
	fours = 0
	for fragid in frag_occurrence[combo]:
		if frag_occurrence[combo][fragid] == 1:
			ones += 1
		elif frag_occurrence[combo][fragid] == 2:
			twos += 1
			sys.stderr.write("%s\t%s\t%s\t%s\t%s\n" % (fragid[0], fragid[1], fragid[2], frag_occurrence[combo][fragid], combo))
		elif frag_occurrence[combo][fragid] == 3:
			threes += 1
			sys.stderr.write("%s\t%s\t%s\t%s\t%s\n" % (fragid[0], fragid[1], fragid[2], frag_occurrence[combo][fragid], combo))
		elif frag_occurrence[combo][fragid] >= 4:
			fours += 1
                        sys.stderr.write("%s\t%s\t%s\t%s\t%s\n" % (fragid[0], fragid[1], fragid[2], frag_occurrence[combo][fragid], combo))			
	print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (frac_human, frac_mouse, bc[combo]["human"], bc[combo]["mouse"], total, cnt[combo], bc1, bc2, "True", "All", ones, twos, threes, fours)
	print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (frac_human_random, frac_mouse_random, bc_random[combo]["human"], bc_random[combo]["mouse"], total_random, cnt[combo], bc1, bc2, "Randomized", "All", ones, twos, threes, fours)
	total = bc_long[combo]["human"] + bc_long[combo]["mouse"]
	if total == 0: continue
	frac_human = float(bc_long[combo]["human"]) / float(total)
	frac_mouse = float(bc_long[combo]["mouse"]) / float(total)
	total_random = bc_random_long[combo]["human"] + bc_random_long[combo]["mouse"]
	if total_random == 0:
		frac_human_random = 0
		frac_mouse_random = 0
	else:
		frac_human_random = float(bc_random_long[combo]["human"]) / float(total_random)
		frac_mouse_random = float(bc_random_long[combo]["mouse"]) / float(total_random)
	print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (frac_human, frac_mouse, bc_long[combo]["human"], bc_long[combo]["mouse"], total, cnt[combo], bc1, bc2, "True", "Long", ones, twos, threes, fours)
	print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (frac_human_random, frac_mouse_random, bc_random_long[combo]["human"], bc_random_long[combo]["mouse"], total_random, cnt[combo], bc1, bc2, "Randomized", "Long", ones, twos, threes, fours)
