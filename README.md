# RCSAnalysisFromT21
Python tool for perform relative critical set analysis on ADF Tape21 files

## To use
Call
	
	python /path/to/RCSAnalysisFromT21.py /path/to/input.txt
	
	
Where the contents of input.txt look like:

	results*.t21	[[1,2],[2,3]]	["path 1","path 2"]
	
so each line is a relative/absolute path to a Tape21 file, a N x 2 list of cp numbers indicating the N paths between cps that should be analyzed, and a N x 1 list of strings to give names to the paths [optional], all separated by single tabs.