import sys
import os



file_list = sys.argv[1]


run = file_list.split('_')[0]

print("run", run)

#outdir = "/pnfs/sbnd/scratch/users/aantonak/PMTAna/"
#outdir = "/pnfs/sbnd/scratch/users/aantonak/PMTAnaTrig/"
outdir = "/pnfs/sbnd/scratch/users/aantonak/PMTAnaOffBeam/"


count = 0
with open(file_list, 'r') as file:
	# loop through each file in the file list
	for line in file:
		new_line = line.strip()	
		#print(new_line)
		output = outdir + "pmt_file"+str(count)+"_"+run+".root"
		count += 1
		#print(output)
		os.system("lar -c run_WaveformScraper.fcl -s "+new_line + " -T "+output)



"""

with open(daq_dir+test_file, 'r') as file:
    # Loop through each line in the file
    for line in file:
        # Strip the newline character and split the line by whitespace
        temp_line = line.strip().split()
        if len(temp_line) > 0:
            if temp_line[0] == "[":
                #print(temp_line)
                i = temp_line.index('170,')
                new_line = temp_line
                new_line[i] = "180,"
                #new_lines.append(new_line)
                final_line = ""
                for w in new_line:
                    final_line += w
                    final_line += " "
                final_line += "\n"
                new_lines.append(final_line)
            else:
                new_lines.append(line)
        else:
            new_lines.append(line)
        #words = line.strip().split()
        # Print the list of words
        print(line)

"""



