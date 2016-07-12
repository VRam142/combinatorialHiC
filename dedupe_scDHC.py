#Goal: quick script to deduplicate sorted bed file
#based on a small window around start and end coordinates
#Input file must be sorted
import os,sys,re

fhi = open(sys.argv[1]) #Open filehandle (sorted BED file)

dedupe = {}
coords_prev = [1,2,3,4,5,6,7,8,9,10] # Initialize a list of previous coordinates
counter = 1
line_prev = ""
for line in fhi:
        split = line.split('\t')
        coords = (split[0], split[1],split[2], split[3],split[4],split[5], split[11], split[12], split[13], split[15]) #Store mapped coordinates and barcode information
        if coords_prev[0] == coords[0] and coords_prev[3] == coords[3] and coords_prev[6] == coords[6] and coords_prev[7] == coords[7]:
#		print coords
#		print coords_prev
                #If chromosomes are the same, check if mapping coordinates fall within 5 bp of each other
    		if int(coords_prev[1]) - 5 <= int(coords[1]) <= int(coords_prev[1]) + 5:
    			if int(coords_prev[2]) - 5 <= int(coords[2]) <= int(coords_prev[2]) + 5:
    				if int(coords_prev[4]) - 5 <= int(coords[4]) <= int(coords_prev[4]) + 5:
    					if int(coords_prev[5]) - 5 <= int(coords[5]) <= int(coords_prev[5]) + 5:
#						print coords_prev
#						print coords_prev
    						coords_prev = coords
    						line_prev = line
    						counter += 1
    						continue
		if coords_prev[8] == coords[8] and coords_prev[9] == coords[9]:
			coords_prev = coords
			line_prev = line
			counter += 1
			continue
		print "%s\t%s" % (line_prev.rstrip('\n'), counter)
		coords_prev = coords
		line_prev = line
		counter += 1
			
        else:
            if coords_prev[0] == 1:
                coords_prev = coords
                line_prev = line
                continue #Ignore print statement if the values are still the same as the initialized list
            print "%s\t%s" % (line_prev.rstrip('\n'), counter)
            coords_prev = coords
            line_prev = line
            counter = 1
