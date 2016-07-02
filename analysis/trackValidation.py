from __future__ import division
import os, sys, argparse, os.path
import ROOT
import array
import math

from pyLCIO import IOIMPL
from pyLCIO import UTIL
from pyLCIO import EVENT
from pyLCIO import IOIMPL

from ROOT import *

import commonFunctions

def update_progress(progress):
	barLength = 20 # Modify this to change the length of the progress bar
	status = ""
	if isinstance(progress, int):
		progress = float(progress)
	if not isinstance(progress, float):
		progress = 0
		status = "error: progress var must be float\r\n"
	if progress < 0:
		progress = 0
		status = "Halt...\r\n"
	if progress >= 1:
		progress = 1
		status = "Done...\r\n"
	block = int(round(barLength*progress))
	text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
	sys.stdout.write(text)
	sys.stdout.flush()

def get_MCParticle_details(MCparticle):
	#MCMom = math.sqrt(math.pow(MCparticle.getMomentum()[0],2) + math.pow(MCparticle.getMomentum()[1],2) + math.pow(MCparticle.getMomentum()[2],2))
	MCMom = math.sqrt(math.pow(MCparticle.getMomentum()[0],2) + math.pow(MCparticle.getMomentum()[1],2))
	MCtheta = 90 - abs(math.degrees(math.atan((MCparticle.getMomentum()[2])/(math.sqrt(math.pow(MCparticle.getMomentum()[0],2) + math.pow(MCparticle.getMomentum()[1],2))))))
	return MCMom, MCtheta

def get_track_details(track):
	trackTransMom = math.fabs((((2.99792*math.pow(10,11))*(math.pow(10,-15))*5)/track.getOmega()))
	trackMom = trackTransMom*math.sqrt(1+pow(track.getTanLambda(),2))
	return trackMom

def track_validation(fileName, numEvents):
	reader = IOIMPL.LCFactory.getInstance().createLCReader()
	reader.open(fileName)	
	valData = []
	count = 0
	badLinkTot = 0
	pionCount = 0
	for event in reader:
		update_progress(count / numEvents)
		count += 1
		MCParticles = event.getCollection("MCParticlesSkimmed")
		truthLink = event.getCollection("TrackMCTruthLink")
		for MCparticle in MCParticles:
			if (MCparticle.getPDG() == 211 or MCparticle.getPDG() == -211):
				pionCount += 1
				MCmom, MCtheta = get_MCParticle_details(MCparticle)
				linkCount = 0
				for link in truthLink:
					if link.getTo() == MCparticle and link.getWeight() == 1.0:
						linkCount += 1
						trackMom = get_track_details(link.getFrom())
						valData.append(['def', MCmom, MCtheta, trackMom - MCmom]) # SHOULD CHANGE THIS SO IT JUST HAS THE WEIGHT AT THE FRONT!!!!!
					if link.getTo() == MCparticle and link.getWeight() != 1.0:
						linkCount += 1
						trackMom = get_track_details(link.getFrom())
						valData.append(['partial', MCmom, MCtheta, trackMom - MCmom])

				if linkCount == 0:
					valData.append(['false', MCmom, MCtheta, 0])
				if linkCount > 1:
					for badLink in range(linkCount):
						badLinkTot += 1
						valData.pop()

	print "----Found " + str(badLinkTot) + " badLinks in file----"
	return valData, pionCount

def maximum_mom(valData):
	maxMom = 0
	for particle in valData:
		if particle[1] > maxMom:
			maxMom = particle[1]
	print "Maximum Momentum = " + str(maxMom)
	return maxMom

def twoDimensionGraphs(name, dataArray, maxMom):
	size = len(dataArray)
	c = TCanvas("c", name)
	graph = ROOT.TGraph2DErrors()
	count = 0
	for mom in range(size):
		for theta in range(size):
			xCoord = mom*(maxMom/size) + ((maxMom/size)/2)
			yCoord = theta*(90/size) + ((90/size)/2)
			graph.SetPoint(count, xCoord, yCoord, dataArray[mom][theta][0])
			graph.SetPointError(count, 0, 0, dataArray[mom][theta][1])
			count += 1
	gStyle.SetPalette(1)
	graph.SetFillColor(29)
	graph.SetMarkerSize(0.8)
	graph.SetMarkerStyle(20)
	graph.SetMarkerColor(kRed)
	graph.SetLineColor(kBlue-3)
	graph.SetLineWidth(2)
	graph.GetXaxis().SetTitle("P_{t}[GeV]")
	graph.GetYaxis().SetTitle("#theta[#circ]")
	#graph.Draw("p0 err")
	graph.Draw("SURF3")
	c.Modified()
	c.Update()
	raw_input()
	img = ROOT.TImage.Create()
	img.FromPad(c)
	img.WriteImage(name + ".png")

def oneDimensionGraphs(name, dataArray, maxMom, cType):
	size = len(dataArray)
	graph = ROOT.TGraphErrors()
	count = 0
	for index in range(size):
		if cType == "P_{t}[GeV]":
			coord = index*(maxMom/size) + ((maxMom/size)/2)
		if cType == "#theta[#circ]":
			coord = index*(90/size) + ((90/size)/2)
		graph.SetPoint(count, coord, dataArray[index][0])
		graph.SetPointError(count, 0, dataArray[index][1])
		count += 1

	if cType == "P_{t}[GeV]":
		graph.SetMinimum(0.75)
		graph.SetMaximum(1)
		graph.GetXaxis().SetLimits(0, 50)
	if cType == "#theta[#circ]":
		graph.SetMinimum(0)
		graph.SetMaximum(1)
		graph.GetXaxis().SetLimits(0, 90)

	graph.GetXaxis().SetTitle(cType)
	graph.GetYaxis().SetTitle("Tracking Efficiency")
	c = TCanvas("C", name)
	c.SetGrid()
	graph.Draw("ALP")
	raw_input()
	img = ROOT.TImage.Create()
	img.FromPad(c)
	img.WriteImage(name + ".png")
	del c, img

def analysisArrayPart(analysisArray, index1, index2):
	defNum, parNum, falNum = 0, 0, 0
	if len(analysisArray[index1][index2]) > 0:
		for part in analysisArray[index1][index2]:
			if part == 'def':
				defNum += 1
			if part == 'partial':
				parNum += 1				
			if part == 'false':
				falNum += 1
			else:
				"Error: This should not be here!!!"
	return defNum, parNum, falNum

def track_analysis(valData, arraySize):
	analysisArray = [[[] for x in range(arraySize)] for x in range(arraySize)]
	#maxMom = maximum_mom(valData)
	maxMom = 40
	momStep = maxMom / arraySize
	angleStep = 90 / arraySize
	#
	limitCount = 0
	for particle in valData:
		if particle[1] < maxMom:
			momIndex = int(math.floor(particle[1]/momStep))
			thetaIndex = int(math.floor(particle[2]/angleStep))
			analysisArray[momIndex][thetaIndex].append(particle[0])
	#
	defArray2d = [[[0, 0] for x in range(arraySize)] for x in range(arraySize)]
	parArray2d = [[[0, 0] for x in range(arraySize)] for x in range(arraySize)]
	falArray2d = [[[0, 0] for x in range(arraySize)] for x in range(arraySize)]
	defMomArray = [[0, 0] for x in range(arraySize)]
	parMomArray = [[0, 0] for x in range(arraySize)]
	falMomArray = [[0, 0] for x in range(arraySize)]
	defThetaArray = [[0, 0] for x in range(arraySize)]
	parThetaArray = [[0, 0] for x in range(arraySize)]
	falThetaArray = [[0, 0] for x in range(arraySize)]

	for index1 in range(arraySize):
		momDef, momPar, momFal, momCount = 0, 0, 0, 0
		thetaDef, thetaPar, thetaFal, thetaCount = 0, 0, 0, 0
		for index2 in range(arraySize):
			defNum, parNum, falNum = analysisArrayPart(analysisArray, index1, index2)

			size = len(analysisArray[index1][index2])
			# print "Def=" + str(defNum / size) + "Par=" + str(parNum / size) + "fal=" + str(falNum / size) + "tot=" + str(size)
			if size > 0:
				# print "Def=" + str(defNum / size) + "Par=" + str(parNum / size) + "fal=" + str(falNum / size) + "tot=" + str(size)
				defArray2d[index1][index2][0] = defNum / size
				parArray2d[index1][index2][0] = parNum / size
				falArray2d[index1][index2][0] = falNum / size
				defArray2d[index1][index2][1] = math.sqrt(((defNum / size)*(1-(defNum / size)))/(size))
				parArray2d[index1][index2][1] = math.sqrt(((parNum / size)*(1-(parNum / size)))/(size))
				falArray2d[index1][index2][1] = math.sqrt(((falNum / size)*(1-(falNum / size)))/(size))

			momDef += defNum
			momPar += parNum
			momFal += falNum

			defNum, parNum, falNum = analysisArrayPart(analysisArray, index2, index1)

			thetaDef += defNum
			thetaPar += parNum
			thetaFal += falNum

			momCount += len(analysisArray[index1][index2])
			thetaCount += len(analysisArray[index2][index1])

		if momCount > 0:
			defMomArray[index1][0] = momDef / momCount
			parMomArray[index1][0] = momPar / momCount
			falMomArray[index1][0] = momFal / momCount
			defMomArray[index1][1] = math.sqrt(((momDef / momCount)*(1-(momDef / momCount)))/(momCount))
			parMomArray[index1][1] = math.sqrt(((momPar / momCount)*(1-(momPar / momCount)))/(momCount))
			falMomArray[index1][1] = math.sqrt(((momFal / momCount)*(1-(momFal / momCount)))/(momCount))
		if thetaCount > 0:
			defThetaArray[index1][0] = thetaDef / thetaCount
			parThetaArray[index1][0] = thetaPar / thetaCount
			falThetaArray[index1][0] = thetaFal / thetaCount
			defThetaArray[index1][1] = math.sqrt(((thetaDef / thetaCount)*(1-(thetaDef / thetaCount)))/(thetaCount))
			parThetaArray[index1][1] = math.sqrt(((thetaPar / thetaCount)*(1-(thetaPar / thetaCount)))/(thetaCount))
			falThetaArray[index1][1] = math.sqrt(((thetaFal / thetaCount)*(1-(thetaFal / thetaCount)))/(thetaCount))

	twoDimensionGraphs("input/validationPixels/partial2D", parArray2d, maxMom)
	twoDimensionGraphs("input/validationPixels/definite2D", defArray2d, maxMom)
	twoDimensionGraphs("input/validationPixels/false2D", falArray2d, maxMom)

	oneDimensionGraphs("input/validationPixels/definite1DMom", defMomArray, maxMom, "P_{t}[GeV]")
	oneDimensionGraphs("input/validationPixels/partial1DMom", parMomArray, maxMom, "P_{t}[GeV]")
	oneDimensionGraphs("input/validationPixels/false1DMom", falMomArray, maxMom, "P_{t}[GeV]")
	oneDimensionGraphs("input/validationPixels/definite1DTheta", defThetaArray, maxMom, "#theta[#circ]")
	oneDimensionGraphs("input/validationPixels/partial1DTheta", parThetaArray, maxMom, "#theta[#circ]")
	oneDimensionGraphs("input/validationPixels/false1DTheta", falThetaArray, maxMom, "#theta[#circ]")

def main():
	# main
	args = commonFunctions.simple_parse_args()

	if not args:
		print 'Invalid Arguments'
		sys.exit(1)

	fileList = commonFunctions.input_files(args.inputDir, ".slcio")
	fileCount = 0
	dataDump = []
	pionCount = 0
	for fileName in fileList:
		fileCount += 1
		print "\nFile -> " + str(fileCount)
		outputName, eventCount, detector, energy = commonFunctions.get_fileData(fileName)
		valData, FilePionCount = track_validation(fileName, eventCount)
		dataDump.extend(valData)
		pionCount += FilePionCount

	print "[Track Validation Done]"
	track_analysis(dataDump, 50)
	print "Total Pions = " + str(pionCount)
	print "[DONE]"

if __name__=='__main__':
	main()