#!/opt/local/bin/python

"""
Need to start Eclipse from command line (using executable inside app) to get
ADF env vars working right.
"""

import kf
import numpy as np
from random import randint

# Takes the path to a Tape21 file WITH BADER CP INFO
# Returns a list of triplets with XYZ coords of each CP
CurrentDir = ''
SubDir = '/Output' + str(randint(1,999999))
PathCoords = []
T21BaseName = ''
ScriptDir = ''

B2A = 0.52917721067

def GetCPCoords(T21FileName):
	T21 = kf.kffile(T21FileName)
	FullList = T21.read("Properties", "CP coordinates")
	NumCPs = T21.read("Properties","CP number of")
	T21.close()
	
	CoordList = []
	
	for i in range(NumCPs):
		tmp = np.empty([3])
		k = 0
		for j in range(i,NumCPs * 3,NumCPs):
			tmp[k] = FullList[j]
			k += 1
		CoordList.append(tmp)
	
	return CoordList
	
# Takes the path to a Tape21 file WITH BADER CP INFO
# Returns a list of CP type numbers
#The "type" numbers go as {1,2,3,4} = {n,c,b,r}
def GetCPTypes(T21FileName):
	TypeList = []
	T21 = kf.kffile(T21FileName)
	Types = T21.stringData("Properties","CP code number for (Rank,Signatu")
	Types = Types.split('\n')
	Types = Types[3:len(Types)-1]
	for i in Types:
		for j in i.split():
			TypeList.append(int(float(j)))
	
	T21.close()
	return TypeList


def GenerateCoords(Pt1, Pt2, NumPts):
	StepDir = Pt2.__sub__(Pt1)
	StepSize = np.linalg.norm(StepDir)
	StepDir /= StepSize
	StepSize /= (NumPts - 1)
	CoordList = [Pt1.__add__(np.multiply(StepDir,Step*StepSize)) for Step in range(NumPts)]
	
	return CoordList


def DensfForCoords(T21FileName, GridCoords):
	import StringIO
	import os
	import subprocess
	global CurrentDir, T21BaseName, ScriptDir
	T41FileName = CurrentDir + '/' + "Tape41_" + str(randint(1,999999)) + ".t41"
	DensfScript = StringIO.StringIO()
	DensfScript.write('''$ADFBIN/densf << eor

INPUTFILE ''' + T21FileName + ''' 
OUTPUTFILE ''' + T41FileName + ''' 

GRID Inline
''')
	
	for i in GridCoords:
		for j in i:
			DensfScript.write(str(j) + ' ')
		DensfScript.write('\n')
	
	DensfScript.write('''End

Density SCF
DenGrad
DenHess

eor''')
	
	DensfScript = DensfScript.getvalue()
# 	
# 	print DensfScript
	
	subprocess.Popen(DensfScript, shell=True, stdout=subprocess.PIPE).stdout.read()
	try:
		os.remove(ScriptDir + '/logfile')
	except: pass
	subprocess.Popen('rm -rf ' + ScriptDir + '/IINFO*' + ScriptDir + '/SINFO*', shell=True, stdout=subprocess.PIPE).stdout.read()
	return T41FileName


def GetDataFromT41(T41FileName):
	T41Data = {}
	T41VarStrs = ['Density']
	for i in ['X','Y','Z']:
		T41VarStrs.append('DensityGrad' + i)
	for i in ['XX','XY','XZ','YY','YZ','ZZ']:
		T41VarStrs.append('DensityHess' + i)
		
	T41 = kf.kffile(T41FileName)
	for i in T41VarStrs:
		T41Data[i] = T41.read('SCF', i)
		
	T41.close()
	T41.delete()
	return T41Data

# HessInd = [
# 	['XX','XY','XZ'],
# 	['XY','YY','YZ'],
# 	['XZ','YZ','ZZ']]
HessInd = [
	['XX','XY','XZ'],
	['XY','YY','YZ'],
	['XZ','YZ','ZZ']]
def GetHessForT41Index(T41Data,Ind):
	H = np.empty([3,3])
	for i in range(3):
		for j in range(3):
			H[i,j] = T41Data['DensityHess' + HessInd[i][j]][Ind]
	return H

GradInd = ['X','Y','Z']
def GetGradForT41Index(T41Data,Ind):
	G = np.empty([3])
	for i in range(3):
		G[i] = T41Data['DensityGrad' + GradInd[i]][Ind]
	return

def GetEigenSystemForT41(T41Data):
	from numpy import linalg as la
	NumPts = len(T41Data['Density'])
	T41Data['EigVals'] = [None]*NumPts
	T41Data['EigVecs'] = [None]*NumPts
	for i in range(NumPts):
		T41Data['EigVals'][i], T41Data['EigVecs'][i] = la.eigh(GetHessForT41Index(T41Data, i))
# 		T41Data['EigVecs'][i] = np.transpose(T41Data['EigVecs'][i])
	return


def tex_escape(text):	
	"""
	    :param text: a plain text message
	    :return: the message escaped to appear correctly in LaTeX
	"""
	import re
	conv = {
	    '&': r'\&',
	    '%': r'\%',
	    '$': r'\$',
	    '#': r'\#',
	    '_': r'\_',
	    '{': r'\{',
	    '}': r'\}',
	    '~': r'\textasciitilde{}',
	    '^': r'\^{}',
	    '\\': r'\textbackslash{}',
	    '<': r'\textless',
	    '>': r'\textgreater',
	}
	regex = re.compile('|'.join(re.escape(unicode(key)) for key in sorted(conv.keys(), key = lambda item: - len(item))))
	return regex.sub(lambda match: conv[match.group()], text)


def RCSMake1DPlot(T41Data, Tol, OutFileName, Res):
	import matplotlib.pyplot as plt
# 	import pylab as plt
# 	from matplotlib import rc
	
	ValList = [[[0],[0]],[[0],[0]]]
	NumPts = len(T41Data['Density'])
	
	
	global CurrentDir, SubDir
	import subprocess
	subprocess.Popen('mkdir -p "' + CurrentDir + SubDir + '"', shell=True, stdout=subprocess.PIPE).stdout.read()
	
	OutCSVName = CurrentDir + SubDir + "/" + OutFileName.replace(' ','_')
	OutCSVName = OutCSVName.rpartition('_(')[0]
	OutCSV = open(OutCSVName + '.csv','w')
	Headings = ["PtNum","X","Y","Z","Rho","Lap"]
	for i in ['1','2','3']:
		Headings.append('DotPdt_' + i)
	for i in ['1','2','3']:
		Headings.append('EigVal_' + i)
	for i in ['X','Y','Z']:
		Headings.append('Grad_' + i)
	for i in ['1','2','3']:
		for j in ['1','2','3']:
			Headings.append('EigVec_' + i + '_' + j)
	for i in ['XX','XY','XZ','YY','YZ','ZZ']:
		Headings.append('Hess_' + i)
			
	for i in Headings:
		OutCSV.write(i + ',')
	OutCSV.write('\n')
	
	
	RCS = [# [[dot products that should be below tolerence],[signs of eigenvalues],[value in plot]]
		[[1,2,3],[1,1,1],[6,6]],
		[[1,2],[-1,-1],[1,2]],
		[[2,3],[1,1],[5,4]],
		[[1,3],[1,1],[3,3]]
		]
	
	for i in range(NumPts):
		G = [T41Data['DensityGrad' + j][i] for j in GradInd]
		E = [T41Data['EigVals'][i][j] for j in range(3)]
# 		EV = np.transpose([T41Data['EigVecs'][i][j] for j in range(3)])
		EV = np.transpose(T41Data['EigVecs'][i])
		DP = [np.absolute(np.dot(G,j)) for j in EV]
		
		PtNum = float(i + 1)
		
		OutCSV.write(str(i+1) + ',')
		for j in [T41Data['Coords'][i]]:
			for k in j:
				OutCSV.write(str(k) + ',')
		OutCSV.write(str(T41Data['Density'][i]) + ',')
		OutCSV.write(str(sum(E)) + ',')
		
		for j in [DP,E,G]:
			for k in j:
				OutCSV.write(str(k) + ',')
		
		for j in EV:
			for k in j:
				OutCSV.write(str(k) + ',')
		for j in [T41Data['DensityHess' + k][i] for k in ['XX','XY','XZ','YY','YZ','ZZ']]:
			OutCSV.write(str(j) + ',')
		OutCSV.write('\n')	
# 		HasVal = False
# 		SameSign = 0
		
		for iRCS in RCS:
			HasVal = True
			for iDP in range(3):
				if ((iDP + 1) in iRCS[0] and DP[iDP] > Tol) or ((iDP + 1) not in iRCS[0] and DP[iDP] <= Tol):
					HasVal = False
					break
			if HasVal:
				RCSInd = 0
				for ei,es in zip(iRCS[0],iRCS[1]):
					if np.sign(E[ei-1]) != es:
						RCSInd = 1
						break
				SignInd = 0 if np.sign(E[iRCS[0][0]-1]) == np.sign(E[iRCS[0][1]-1]) else 1
				ValList[SignInd][0].append(PtNum) 
				ValList[SignInd][1].append(iRCS[2][RCSInd])
				
				
				
		
# 		# CP (not used)
# 		if DP[0] <= Tol and DP[1] <= Tol and DP[2] <= Tol:
# 			Val = 6
# 			HasVal = True
# 			
# 		# Ridge and r-connector
# 		elif DP[0] <= Tol and DP[1] <= Tol and DP[2] > Tol:
# 			HasVal = True
# 			if E[0] < 0:
# 				if E[1] < 0: Val = 1
# 				else:
# 					SameSign = 1
# 					Val = 2
# 			else:
# 				Val = 2
# 				if E[1] < 0: SameSign = 1
# 					
# 		# Valley and v-connector
# 		elif DP[1] <= Tol and DP[2] <= Tol and DP[0] > Tol:
# 			HasVal = True
# 			if E[1] > 0:
# 				if E[2] > 0: Val = 5
# 				else:
# 					SameSign = 1
# 					Val = 4
# 			else:
# 				Val = 4
# 				if E[2] > 0: SameSign = 1
# 					
# 		# m-connector
# 		elif DP[0] <= Tol and DP[2] <= Tol and DP[1] > Tol:
# 			HasVal = True
# 			Val = 3
# 			if np.sign(E[2]) != np.sign(E[0]): SameSign = 1 
# 			
# 		if HasVal:
# 			ValList[SameSign][0].append(PtNum)
# 			ValList[SameSign][1].append(Val)
	
	OutCSV.close()
				
	for j in ValList: 
		j[0].append(NumPts)
		j[1].append(0)
		
	NumPts = float(NumPts)
	
	ValList = [[[j / NumPts * 100. for j in k[0]],[j for j in k[1]]] for k in ValList]
	
	xVals = range(0,125,25)
	yVals = range(1,7)
	
	fig = plt.figure()
	ax = fig.add_subplot(111) # subplot required to change the aspect ratio
	ax.set_aspect(10)
	ax.plot(ValList[0][0], ValList[0][1], 'bs', ValList[1][0], ValList[1][1], 'r+')
	
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.rc('font', size=16)
	yLabels = ['Ridge',r'$R$-connector',r'$M$-connector',r'$V$-connector','Valley','Critical point']
	
	plt.xticks(xVals, xVals)
	plt.yticks(yVals,yLabels,rotation=0)
	plt.axis([-2,103,0.5,5.5])
	plt.xlabel('\% of path langth')
	
	fig.tight_layout()
	
	fig.savefig(CurrentDir + SubDir + "/" + OutFileName.replace(' ','_') + '.png', dpi=Res)
	plt.axis([-2,103,0.5,6.5])
	fig.savefig(CurrentDir + SubDir + "/" + OutFileName.replace(' ','_') + '_wCP.png', dpi=Res)
	
	plt.title(tex_escape(OutFileName.replace('_',' ')))
	plt.axis([-2,103,0.5,5.5])
	fig.savefig(CurrentDir + SubDir + "/" + OutFileName.replace(' ','_') + '_wTitle.png', dpi=Res)
	plt.axis([-2,103,0.5,6.5])
	fig.savefig(CurrentDir + SubDir + "/" + OutFileName.replace(' ','_') + '_wTitle_wCP.png', dpi=Res)
	
# 	plt.show()
	
	plt.close()
					
	return

def GetCPsFromCoords(CPList,NewCoords,NewTypes,OldCoords,OldTypes):
	CPNums = []
	
	for CP in CPList:
		MinDist = 1e100
		MinInd = -1
		for i in range(len(NewCoords)):
			if NewTypes[i] == OldTypes[CP-1]:
				TmpDist = np.linalg.norm(np.subtract(NewCoords[i],OldCoords[CP-1]))
				if TmpDist < MinDist:
					MinDist = TmpDist
					MinInd = i
		if MinInd >= 0:
			CPNums.append(MinInd + 1)
	
	return CPNums
	
	
CPTypeInd = ['N','C','B','R']
def main():
	import ast
	import os
	import sys
	import subprocess
	global T21BaseName, ScriptDir, CurrentDir
	
	Res = 300 # dpi of resulting plots
	
	NumPts = 301
	Tols = [10**(-e) for e in range(3,9)]
# 	InputFilePath = "/Users/Haiiro/Safe/tmp/input.txt"
# 	InputFilePath = "/Users/Haiiro/ADFdata/1RidgePaper/first_round_data/input.txt"
# 	InputFilePath = "/Users/Haiiro/ADFdata/1RidgePaper/first_round_data/input_copy.txt"
	InputFilePath = "/Volumes/Kuroi1_Misc/workspace/1RidgePaper/first_round_data/input.txt"
	
	if len(sys.argv) > 1:
		InputFilePath = sys.argv[1]
	if len(sys.argv) > 2:
		NumPts = int(sys.argv[2])
	if len(sys.argv) > 3:
		Tols = sys.argv[3]
		if len(Tols.split(',')) > 1: 
			Tols = Tols.split(',')
			Tols = [10**(-e) for e in range(int(Tols[0]),int(Tols[1])+1)]
		else: Tols = [10**(-int(Tols))]
		

	InputFile = open(InputFilePath, 'r')
	InputFileContents = InputFile.readlines()
	InputFile.close()
	
	ScriptDir = os.path.dirname(os.path.realpath(__file__))
	InputDir = os.path.dirname(os.path.realpath(InputFilePath))
	
	FileNum = 0
	
	for aLine in InputFileContents:
		FileNum += 1
		print '\nRunning file ' + str(FileNum) + ' of ' + str(len(InputFileContents))
		
		aLine = aLine.split('\t')
		if len(aLine) > 0:
			InputT21FileName = aLine[0]
			TargetCPs = ast.literal_eval(aLine[1])
			
			print "File: %s with CPs: %s" % (InputT21FileName, str(TargetCPs))
			
			T21BaseName = os.path.basename(InputT21FileName)
				
			if T21BaseName != InputT21FileName:
				CurrentDir = os.path.dirname(InputT21FileName)
			else:
				CurrentDir = InputDir
				
			InputT21FileName = CurrentDir + '/' + T21BaseName
			
			print "Adjusted file name: %s" % InputT21FileName
			
			if '?' in InputT21FileName or '*' in InputT21FileName:
				T21FileNames = subprocess.Popen("ls " + InputT21FileName, shell=True, stdout=subprocess.PIPE).stdout.read()
				T21FileNames = T21FileNames.split('\n')
			else:
				T21FileNames = [InputT21FileName]
				
			print "All T21s for file: %s" % str(T21FileNames)
				
			OrigCPCoords = []
			OrigCPTypes = []
			
			for T21FileName in T21FileNames:
				print "Running T21: %s" % T21FileName
				if os.path.exists(T21FileName):
					T21BaseName = os.path.basename(T21FileName)
					
					if '.' in T21BaseName:
						T21BaseName = T21BaseName.rpartition('.')[0]
						
# 					T21FileName = "'" + T21FileName + "'"
				
					CPCoords = GetCPCoords(T21FileName)
					CPTypes = GetCPTypes(T21FileName)
					
					if len(OrigCPCoords) == 0: OrigCPCoords = CPCoords
					if len(OrigCPTypes) == 0: OrigCPTypes = CPTypes
					
					OutFileNames = []
					
					if len(aLine) > 2:
						OutFileNames = ast.literal_eval(aLine[2])
						if len(OutFileNames) != len(TargetCPs):
							OutFileNames = []
							break
						Tmp = OutFileNames
						OutFileNames = []
						for i in Tmp:
							TmpStr = T21BaseName + '_' + i + '_'
# 							First = True
# 							for i in c:
# 								if First: First = False
# 								else: TmpStr += '-'
# 								TmpStr += CPTypeInd[CPTypes[i-1]-1] + str(i)
							OutFileNames.append(TmpStr)
					
					if len(OutFileNames) == 0:
						for c in TargetCPs:
							TmpStr = T21BaseName + '_'
# 							First = True
# 							for i in c:
# 								if First: First = False
# 								else: TmpStr += '-'
# 								TmpStr += CPTypeInd[CPTypes[i-1]-1] + str(i)
							OutFileNames.append(TmpStr)
					
					global PathCoords
					
					for Pts, PathName in zip(TargetCPs, OutFileNames):
						CPNums = GetCPsFromCoords(Pts, CPCoords, CPTypes, OrigCPCoords, OrigCPTypes)
						
						print "CP pair %s (%s)" % (str(CPNums), Pts)
						
						OutFileName = PathName
						
						First = True
						for i in CPNums:
							if First: First = False
							else: OutFileName += '-'
							OutFileName += CPTypeInd[CPTypes[i-1]-1] + str(i)
							
						print "Outfile base: %s" % OutFileName
						
						PathCoords = GenerateCoords(CPCoords[CPNums[0]-1] * B2A, CPCoords[CPNums[1]-1] * B2A, NumPts)
						T41FileName = DensfForCoords(T21FileName, PathCoords)
						T41Data = GetDataFromT41(T41FileName)
						T41Data['Coords'] = PathCoords
						GetEigenSystemForT41(T41Data)
						for Tol in Tols:
							RCSMake1DPlot(T41Data, Tol, OutFileName + '_(' + '%.0e' % Tol + ')', Res)
				# 		break
				else:
					print "File not found: %s" % T21FileName
		
		
	print "Finished"
	return


main()