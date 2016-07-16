# RCSAnalysisFromT21
Python tool for perform relative critical set analysis on ADF Tape21 files

## To use
Call
	
	python /path/to/RCSAnalysisFromT21.py /path/to/input.txt
	
	
Where the contents of input.txt look like:

	results*.t21	[[1,2],[2,3]]	["path 1","path 2"] ["img1.png","img2.png"]
	
so each line is a relative/absolute path to a Tape21 file, a N x 2 list of cp numbers indicating the N paths between cps that should be analyzed, a N x 1 list of strings to give names to the paths [optional], and a N x 1 list of paths of images to go with each CP pair that will be placed to the right of the resulting plots (this uses eval() so you can use a list comprehension if you want, but beware of malicious or faulty code in that eval() statement), all separated by single tabs.