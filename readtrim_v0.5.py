#!usr/bin/env python3
import sys, os, gzip, argparse

################################################################
#### @Author:	Josephine Strange							####
#### @Version:	0.5											####
#### @Date:		09-Dec-2019									####
####														####
####				Read trimmer for NGS data				####
####					  The basics						####
####														####
################################################################
#### Define functions ####
def detect_Q_encoding(baseQualities):
	"""Returns Phred encoding for first unique char encountered (max Q40)"""
	for char in baseQualities:
		# Char only in phred 33 convention
		if ord(char) >= 33 and ord(char) < 64:
		 	return "Phred33"
		 # Char only in phred 64 convention
		elif ord(char) >= 74 and ord(char) <= 104 :
			return "Phred64"
		# If char is ambiguous
		elif ord(char) >= 64 and ord(char) < 74:
			pass
		# Char in neither convention
		else:
			sys.stderr.write("\nUnknown encoding encountered ({} with ASCII value {}). Terminating!\n".format(char, ord(char)))
			sys.exit(1)

def average_Q_score(window):
	""" Returns the Q-score for a given quality window"""
	sumScore = 0
	for i in range(len(window)):
		sumScore += ord(window[i])

	return sumScore/len(window)

################################################################
#### Setup argument parser #### 
# Setup parser object
parser = argparse.ArgumentParser()

# The file parser accepts a list of values - requires 1 or more arguments (nargs = "+")
parser.add_argument("-f", "--file", nargs="+", help="specify file(s): Maximum two of paired-end reads", required=True) # Required optional arg

# Trim ends based on number of nucleotides . requires 0 or 1 arguments (nargs = "?")
parser.add_argument("-t5", "--trim_5", nargs="?", help="number of nucleotides to trim from 5'", type=int)
parser.add_argument("-t3", "--trim_3", nargs="?", help="number of nucleotides to trim from 3'", type=int) 

# Trim based on quality - requires 0 or 1 arguments (nargs = "?")
parser.add_argument("-tq", "--quality",  nargs="?", help="trim reads based on minimum Phred quality", type=int)

# Filter based on average 
parser.add_argument("-fa", "--average", help="minimum average Q-score to retain read after trimming", type=int)

# Type of Phred encoding - rTrue or False (Default False)
parser.add_argument("-phred33", action="store_true", help="phred score encoding 33")
parser.add_argument("-phred64", action="store_true", help="phred score encoding 64")

# Maximum number of unknown nucletides before read is discarded - requires 0 or 1 arguments (nargs = "?")
parser.add_argument("-unk", "--unknown", nargs="?", help="maximum number of unknown (N) bases to retain read", type=int)
parser.add_argument("-l", "--length", nargs=1, help="minimum number of nucleotides to retain read after trimming", type=int)

# Define argument objects
args = parser.parse_args()

################################################################
#### File and argument parsing ####
fileExt = (".fq", ".fastq", ".gz")

# Iterate over list of files from argparse discard files of wrong file format
readFiles = []
for file in args.file: 
	if not file.endswith(fileExt):
		sys.stdout.write("Detected incompatible file format: \"{}\"... Discarding file!\n".format(os.path.splitext(file)[1]))
	else:
		readFiles.append(file)

################################################################
#### Output beautification #### 
#### N threshold ####
if args.unknown is not None:
	nThreshold = args.unknown
elif args.unknown is None:
	nThreshold = "N/A"

#### Trim bases from 5' and/or 3' end ####
if args.trim_5 is not None:
	trim_5 = args.trim_5
elif args.trim_5 is None:
	trim_5 = "N/A"

if args.trim_3 is not None:
	trim_3 = args.trim_3
elif args.trim_3 is None:
	trim_3 = "N/A"

#### Trim reads based on quality ####
trimQuality = False
if args.quality is not None:
	trimQuality = True
	qThreshold = args.quality
elif args.quality is None:
	trimQuality = False
	qThreshold = "N/A"

#### Filter reads on length ####
if args.length is not None:
	minLength = args.length
elif args.length is None:
	minLength = "N/A"

################################################################
#### Work with gzipped/fastq files ####
# Initialize variables
fileCompression, fileWrite = None, True


for file in readFiles:
	removedReads,trimmedReads, trimmedBases = 0, 0, 0
	phredEncoding = None
	sumUnk = 0

	#### Open input read files - open and setup output FASTQ files and log ####
	try:
		# Open and setup log
		log_name = os.path.basename(file).split(".")[0] + ".log"
		log = open(log_name, "w")
		log.write("")
		
		# Output FASTQ filename
		fastq_out_name = os.path.basename(file).split(".")[0] + ".trim.fastq"
		
		# Open and setup compressed FASTQ file
		if file.endswith(".gz"):
			fastq_in = gzip.open(file, "rt")
			sys.stdout.write("\nDetected compressed FASTQ file...")

			# Outfile specifications
			fastq_out_name += ".gz"
			fastq_out = gzip.open(fastq_out_name, "wb")
			fastq_out.write("".encode()) # Override file
			fileCompression = True

		# Open and setup uncompressed FASTQ file
		else:
			fastq_in = open(file, "r")
			sys.stdout.write("\nDetected raw FASTQ file...")

			# Outfile specifications
			fastq_out = open(fastq_out_name, "w")
			fastq_out.write("") # Override file
			fileCompression = False

	except IOError as err:
		sys.stderr.write("\nError when loading file: {}\n".format(str(err)))
		sys.exit(1)
	
	# Define variables
	lineNum = 0
	numReads = 0
	numUnk = 0
	header,sequence,comment,quality = "","","",""

	# Iterate over lines in the current file
	for line in fastq_in:
		lineNum += 1

		################################################################
		# Stores header (1st line of read)
		if lineNum == 1:
			header = line
						
		################################################################
		# Stores read sequence (2nd line of a read)
		elif lineNum == 2:
			sequence = line.rstrip()

			#### if args.unknown has been specified ####
			unkNum = 0
			if args.unknown is not None:
				for nucleotide in sequence:
					if nucleotide is "N":
						unkNum += 1
				sumUnk += unkNum
				# Discard read if number of unknown bases above threshold
				if unkNum > args.unknown:
					fileWrite = False
					removedReads += 1
					
		################################################################
		# Stores comment field (3rd line of read)
		elif lineNum == 3:
			comment = line

		################################################################
		# Stores base qualities (4th line of a read) - entire read (4 lines) in memory
		elif lineNum == 4:
			quality = line.rstrip()
			numReads += 1

			#### Encoding ####
			if phredEncoding is None:
				# If encoding has been specified by user set encoding to user input
				if args.phred33:
					phredEncoding = "Phred33"
				elif args.phred64:
					phredEncoding = "Phred64"

				# If encoding has not been specified in user input - auto-detect
				else:
					sys.stdout.write("\nEncoding not specified by user. Auto-detecting...")
					phredEncoding = detect_Q_encoding(quality)
					sys.stdout.write("\nEncoding detected: {}\n".format(phredEncoding))
			
			#### Trim X nucleotides from either or both ends ####
			if args.trim_3 is not None or args.trim_5 is not None:
				if fileWrite:
					trimmedReads += 1
				
					if args.trim_3 is not None:
						trimmedBases += args.trim_3
						sequence = sequence[:-args.trim_3]
						quality = quality[:-args.trim_3]
	
					if args.trim_5 is not None:
						trimmedBases += args.trim_5
						sequence = sequence[args.trim_5:]
						quality = quality[args.trim_5:]
			
			#### Trim the reads based on average Q-score of moving window #### 
			if args.quality is not None:
				if fileWrite: 
					pass

			#### Filter out reads shorter than user specificed minimum length ####
			if args.length is not None:
				if fileWrite:
					if len(sequence) < args.length:
						fileWrite = False
						removedReads += 1
						
			#### Filter out reads with a mean quality lower than specified after trimming ####
			if args.average is not None:
				if fileWrite:

					sumScore = 0
					for i in range(len(sequence)):
						if phredEncoding == "Phred33":
							sumScore += ord(sequence[i])-33
						elif phredEncoding == "Phred64":
							sumScore += ord(sequence[i])-64
					
					avgScore = sumScore // len(sequence)
					if avgScore < args.average:
						fileWrite = False
						removedReads += 1


			#### Write to fastq file ####
			if fileWrite:
				if fileCompression: 
					fastq_out.write("{}{}\n{}{}\n".format(header,sequence,comment,quality).encode("ascii"))

				else:
					fastq_out.write("{}{}\n{}{}\n".format(header,sequence,comment,quality))
				
			# Re-initialize controller variables 
			lineNum = 0
			fileWrite = True	

		

	#### Write to log ####	
	sys.stdout.write("\nWriting to {}\n".format(log_name))
	
	log.write("\
Input:\t{}\n\
Output:\t{}\n\n\
Reads in input:\t{}\n\
Q-score encoding:\t{}\n\
Q-score threshold:\t{}\n\
Bases trimmed from 5':\t{}\n\
Bases trimmed from 3':\t{}\n\
Trimmed reads (total):\t{}\n\
Trimmed bases (total):\t{}\n\
N threshold:\t{}\n\
Avg. Q-score threshold:\t{}\n\
Removed reads:\t{}\n\
Reads in output:\t{}\n\
	".format(\
		os.path.basename(file),\
		fastq_out_name,\
	 	numReads,\
	 	phredEncoding,\
	 	qThreshold,\
	 	trim_5,\
	 	trim_3,\
	 	trimmedReads,\
	 	trimmedBases,\
	 	nThreshold,\
	 	args.average,\
	 	removedReads,\
	 	(numReads-removedReads)))
		
fastq_in.close()
fastq_out.close()
log.close()