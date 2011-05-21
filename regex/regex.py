import re, sys

# Get inFile name
inFileName = sys.argv[1]

# Open inFile for reading
inFile = open(inFileName, "r")


# Read in the inFile contents as a string
inFileAsStr = inFile.read()


# Close inFile
inFile.close()


# Define & compile the regex patterns
patterns = [r"public:"]
for pattern in patterns:
	pattern = re.compile(pattern, re.S)


# Find the patterns in the inFile string
classIter = re.finditer(patterns[0], inFileAsStr)


# Gets length of classIter
classIterLen = 0
for iter in classIter:
	classIterLen = classIterLen + 1
classIter = classIter._iter_()

# Print number of matches and the matches themselves
print classIterLen, "matches found..."
for match in classIter.last():
	print "================\r"
	print match.group(0)



