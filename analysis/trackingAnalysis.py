from __future__ import division
import os
import sys
import argparse
import os.path
import ROOT
import array
import math

from pyLCIO import IOIMPL
from pyLCIO import UTIL
from pyLCIO import EVENT
from pyLCIO import IOIMPL

from ROOT import *

import commonFunctions

def track_validation(fileName, numEvents):
	reader = IOIMPL.LCFactory.getInstance().createLCReader()
	reader.open(fileName)	
	pionTrackData = []
	count, pionCount = 0, 0
	for event in reader:
		commonFunctions.update_progress(count/numEvents)
		count += 1
		MCParticles = event.getCollection("MCParticlesSkimmed")
		truthLink = event.getCollection("TrackMCTruthLink")
		if truthLink.getNumberOfElements() != event.getCollection("Tracks").getNumberOfElements():
			print "Error: Different number of track elements to truthLink elements!!!"
		for MCparticle in MCParticles:
			if (MCparticle.getPDG() == 211 or MCparticle.getPDG() == -211):
				pionCount += 1
				MCMomTot, MCMomTrans, MCtheta = get_MCParticle_details(MCparticle)
				linkCount = 0
				for link in truthLink:
					if link.getTo() == MCparticle:
						linkCount +=1
						trackTransMom = get_track_details(link.getFrom())
						pionTrackData.append([link.getWeight(), MCMomTot, MCtheta, MCMomTrans, trackTransMom - MCMomTrans])
				if linkCount == 0:
					pionTrackData.append([0, MCMomTot, MCtheta, 0])
				elif linkCount > 1:
					temp = []
					for badLink in range(linkCount):
						temp.append(pionTrackData.pop())
						temp[badLink][0] = linkCount
					pionTrackData.extend(temp)

	return pionTrackData, pionCount

def get_MCParticle_details(MCparticle):
	MCMomTot = math.sqrt(math.pow(MCparticle.getMomentum()[0],2) + math.pow(MCparticle.getMomentum()[1],2) + math.pow(MCparticle.getMomentum()[2],2))
	MCMomTrans = math.sqrt(math.pow(MCparticle.getMomentum()[0],2) + math.pow(MCparticle.getMomentum()[1],2))
	MCtheta = 90 - abs(math.degrees(math.atan((MCparticle.getMomentum()[2])/(math.sqrt(math.pow(MCparticle.getMomentum()[0],2) + math.pow(MCparticle.getMomentum()[1],2))))))
	return MCMomTot, MCMomTrans, MCtheta

def get_track_details(track):
	trackTransMom = math.fabs((((2.99792*math.pow(10,11))*(math.pow(10,-15))*5)/track.getOmega()))
	#trackMom = trackTransMom*math.sqrt(1+pow(track.getTanLambda(),2))
	return trackTransMom

def analysisArrayBin(arrayBin):
	defNum, parNum, falNum, size = 0, 0, 0, 0
	size = len(arrayBin)
	if size > 0:
		for part in arrayBin:
			if part == 1:
				defNum += 1
			if part < 1 and part > 0:
				parNum += 1				
			if part > 1:
				falNum += 1
	return defNum, parNum, falNum, size

def mom_graph(delta_mom, bins, width):
	hist = ROOT.TH1F("momDiff", "momDiff", bins, -width, width)
	for mom in delta_mom:
		hist.Fill(mom)
	hist.Fit("gaus")
	sigma = 0
	sigmaError = 0
	try:
		sigma = hist.GetFunction("gaus").GetParameter(2)
		sigmaError = hist.GetFunction("gaus").GetParError(2)
	except: print "Boom"
	if sigma != 0:
		return sigma, sigmaError
	else:
		return sigma, sigmaError

def analysis_arrays(pionTrackData, bins, maxMom):
	analysisArray = [[[] for x in range(bins)] for x in range(bins)]
	momDiffArray = [[] for x in range(bins)]
	momStep = maxMom / bins
	angleStep = 90 / bins
	for particle in pionTrackData:
		if particle[1] < maxMom:
			momIndex = int(math.floor(particle[1]/momStep))
			thetaIndex = int(math.floor(particle[2]/angleStep))
			analysisArray[momIndex][thetaIndex].append(particle[0])	
			momTransIndex = int(math.floor(particle[3]/momStep))
			if particle[0] == 1:
				momDiffArray[momTransIndex].append(particle[4])	

	efficiencyArray = [[[0, 0, 0, 0] for x in range(bins)] for x in range(bins)] #defNum, parNum, falNum, size
	momSigmaArray = [[0, 0] for x in range(bins)] #Sigma, SigmaError

	for index1 in range(bins):
		sigma, sigmaError = mom_graph(momDiffArray[index1],bins,(index1+1)*0.1)
		momSigmaArray[index1][0] = sigma / pow(((index1*momStep)+(momStep/2)),2)
		momSigmaArray[index1][1] = sigmaError # NEED TO THINK ABOUT THIS ERROR, need to include uncertainty on the transverse momentum.
		for index2 in range(bins):
			defNum, parNum, falNum, size = analysisArrayBin(analysisArray[index1][index2])
			efficiencyArray[index1][index2][0] = defNum
			efficiencyArray[index1][index2][1] = parNum
			efficiencyArray[index1][index2][2] = falNum
			efficiencyArray[index1][index2][3] = size

	return efficiencyArray, momSigmaArray

def graphs_mom(name, efficiencyArray, bins, maxMom):
	thetaCut = math.floor(bins/2)
	BCGraph = ROOT.TGraphErrors()
	ACGraph = ROOT.TGraphErrors()
	for mom in range(bins):
		defNumBellow, defNumAbove, countBellow, countAbove = 0,0,0,0
		for theta in range(bins):
			if theta <= thetaCut:
				defNumBellow += efficiencyArray[mom][theta][0]
				countBellow += efficiencyArray[mom][theta][3]
			if theta > thetaCut:
				defNumAbove += efficiencyArray[mom][theta][0]
				countAbove += efficiencyArray[mom][theta][3]

		coord = mom*(maxMom/bins) + ((maxMom/bins)/2)
		if countBellow > 0:
			BCGraph.SetPoint(mom, coord, defNumBellow / countBellow)
			BCGraph.SetPointError(mom, 0, math.sqrt(((defNumBellow / countBellow)*(1-(defNumBellow / countBellow)))/(countBellow))) 
		else:
			BCGraph.SetPoint(mom, coord, 0)
			BCGraph.SetPointError(mom, 0, 0)
		if countAbove > 0:
			ACGraph.SetPoint(mom, coord, defNumAbove / countAbove)
			ACGraph.SetPointError(mom, 0, math.sqrt(((defNumAbove / countAbove)*(1-(defNumAbove / countAbove)))/(countAbove))) 
		else:
			ACGraph.SetPoint(mom, coord, 0)
			ACGraph.SetPointError(mom, 0, 0) 

	BCGraph.GetXaxis().SetTitle("P[GeV]")
	BCGraph.GetYaxis().SetTitle("Tracking Efficiency")
	BCGraph.SetName(name + "_<45Degrees")
	BCGraph.SetMinimum(0.8)
	BCGraph.SetMaximum(1)
	BCGraph.GetXaxis().SetLimits(0, 50)
	ACGraph.SetName(name + "_>45Degrees")
	return BCGraph, ACGraph

def graphs_theta(name, efficiencyArray, bins, maxMom):
	momCut = 1
	BCGraph = ROOT.TGraphErrors()
	ACGraph = ROOT.TGraphErrors()
	for theta in range(bins):
		defNumBellow, defNumAbove, countBellow, countAbove = 0,0,0,0
		for mom in range(bins):
			if mom <= momCut:
				defNumBellow += efficiencyArray[mom][theta][0]
				countBellow += efficiencyArray[mom][theta][3]

			if mom > momCut:
				defNumAbove += efficiencyArray[mom][theta][0]
				countAbove += efficiencyArray[mom][theta][3]

		coord = theta*(90/bins) + ((90/bins)/2)
		if countBellow > 0:
			BCGraph.SetPoint(theta, coord, defNumBellow / countBellow)
			BCGraph.SetPointError(theta, 0, math.sqrt(((defNumBellow / countBellow)*(1-(defNumBellow / countBellow)))/(countBellow))) 
		else:		
			BCGraph.SetPoint(theta, coord, 0)
			BCGraph.SetPointError(theta, 0, 0) 		
		if countAbove > 0:
			ACGraph.SetPoint(theta, coord, defNumAbove / countAbove)
			ACGraph.SetPointError(theta, 0, math.sqrt(((defNumAbove / countAbove)*(1-(defNumAbove / countAbove)))/(countAbove))) 
		else:
			ACGraph.SetPoint(theta, coord, 0)
			ACGraph.SetPointError(theta, 0, 0)  

	BCGraph.GetXaxis().SetTitle("#theta[#circ]")
	BCGraph.GetYaxis().SetTitle("Tracking Efficiency")
	BCGraph.SetName(name + "_<1Gev")
	BCGraph.SetMinimum(0.5)
	BCGraph.SetMaximum(1)
	BCGraph.GetXaxis().SetLimits(0, 90)	
	ACGraph.SetName(name + "_> 10Gev")

	return BCGraph, ACGraph

def res_graph(inputDir, momSigmaArrayStrips, momSigmaArrayPixels, size, maxMom):
	graphStrips = ROOT.TGraph()
	graphPixels = ROOT.TGraph()
	count = 0
	for index in range(size):
		coord = index*(maxMom/size) + ((maxMom/size)/2)
		graphStrips.SetPoint(count, coord, momSigmaArrayStrips[index][0])
		graphPixels.SetPoint(count, coord, momSigmaArrayPixels[index][0])
		#graph.SetPointError(count, 0, momSigmaArray[index][1])
		count += 1
	c = TCanvas("C", "res")
	c.SetLogy()
	c.SetLogx()
	c.SetGrid()
	graphStrips.GetXaxis().SetTitle("P_{t}[GeV]")
	graphStrips.GetYaxis().SetTitle("#sigma(P_{t}) / P_{t}^{2}")
	graphStrips.SetMinimum(pow(10,-5))
	graphStrips.SetMaximum(pow(10,-2))
	graphStrips.GetXaxis().SetLimits(0, 50)
	graphStrips.SetName("strips")
	graphPixels.SetName("pixels")
	graphStrips.Draw("ALP")
	graphPixels.Draw("SAME")
	c.Modified()
	c.Update()
	raw_input()
	img = ROOT.TImage.Create()
	img.FromPad(c)
	img.WriteImage(inputDir + "res.png")
	del c, img

def parse_args():
	# Takes in arguments from the command line when script is called and returns them to main().
	currentDir = os.getcwd()
	parser = argparse.ArgumentParser(description='Processes .slcio file...'
						  ,epilog='In case of questions or problems, contact jt12194@my.bristol.ac.uk')

	parser.add_argument('-p', '--pixels', help='Input .slcio fileDirectory for pixels')

	parser.add_argument('-s', '--strips', help='Input .slcio fileDirectory for normal strips')

	parser.add_argument('-o', '--outputDirectory',
						help='Output directory for plots to be saved to, default = current directory.',
						default=currentDir)

	parser.add_argument('-bins', '--numberOfBins',
						help='The number of bins for the histograms.',
						default=100)

	return parser.parse_args()

def graphs_2d(name, efficiencyArray, bins, maxMom, inputDir):
	# Takes in an efficiencyArray and produces a 2D efficiency plot from it.
	graph = ROOT.TGraph2DErrors()
	count = 0
	for mom in range(bins):
		for theta in range(bins):
			xCoord = mom*(maxMom/bins) + ((maxMom/bins)/2)
			yCoord = theta*(90/bins) + ((90/bins)/2)
			true = efficiencyArray[mom][theta][0]
			num = efficiencyArray[mom][theta][3]
			if num > 0:
				graph.SetPoint(count, xCoord, yCoord, true / num)
				graph.SetPointError(count, 0, 0, math.sqrt(((true / num)*(1-(true / num)))/(num)))
			else:
				graph.SetPoint(count, xCoord, yCoord, 0)
				graph.SetPointError(count, 0, 0, 0)				
			count += 1

	gStyle.SetPalette(1)
	graph.GetXaxis().SetTitle("P[GeV]")
	graph.GetYaxis().SetTitle("#theta[#circ]")
	graph.SetName(name + "2D")
	c = TCanvas("c", "canvas")
	graph.Draw("SURF3")
	c.Modified()
	c.Update()
	raw_input()
	img = ROOT.TImage.Create()
	img.FromPad(c)
	img.WriteImage(name + "_2D.png")

def plot_1D_graph(name, graphs, outputDirectory):
	c = TCanvas("C", name)
	c.SetGrid()
	graphs[0].Draw("ALP")
	graphs.pop(0)
	for graph in graphs:
		graph.Draw("SAME")

	c.Modified()
	c.Update()
	raw_input()
	img = ROOT.TImage.Create()
	img.FromPad(c)
	img.WriteImage(outputDirectory + name + ".png")

def main():
	args = parse_args()
	if not args:
		print 'Invalid Arguments'
		sys.exit(1)

	bins = 50
	maxMom = 50

	stripsFileList = commonFunctions.input_files(args.strips, ".slcio")
	pixelsFileList = commonFunctions.input_files(args.pixels, ".slcio")

	fileCount, pionCount = 0, 0
	pionTrackDataStrips = []
	print "\n[LOADING STRIP FILES]"
	for fileName in stripsFileList:
		fileCount += 1
		print "\nFile -> " + str(fileCount)
		outputName, eventCount, detector, energy = commonFunctions.get_fileData(fileName)
		filePionTrackData, FilePionCount = track_validation(fileName, eventCount)
		pionTrackDataStrips.extend(filePionTrackData)
		pionCount += FilePionCount
	fileCount, pionCount = 0, 0
	pionTrackDataPixels = []
	print "\n[LOADING PIXEL FILES]"
	for fileName in pixelsFileList:
		fileCount += 1
		print "\nFile -> " + str(fileCount)
		outputName, eventCount, detector, energy = commonFunctions.get_fileData(fileName)
		filePionTrackData, FilePionCount = track_validation(fileName, eventCount)
		pionTrackDataPixels.extend(filePionTrackData)
		pionCount += FilePionCount
	print "\n[LOADING FILES DONW] Total Pions used = " + str(pionCount)

	print "\n[CREATING ARRAYS]"
	efficiencyArrayStrips, momSigmaArrayStrips = analysis_arrays(pionTrackDataStrips, bins, maxMom)
	efficiencyArrayPixels, momSigmaArrayPixels = analysis_arrays(pionTrackDataPixels, bins, maxMom)

	print "\n[CREATING GRAPHS]"
	BCGraphMomStrips, ACGraphMomStrips = graphs_mom("strips", efficiencyArrayStrips, bins, maxMom)
	BCGraphMomPixels, ACGraphMomPixels = graphs_mom("pixels", efficiencyArrayPixels, bins, maxMom)
	BCGraphThetaStrips, ACGraphThetaStrips = graphs_theta("strips", efficiencyArrayStrips, bins, maxMom)
	BCGraphThetaPixels, ACGraphThetaPixels = graphs_theta("pixels", efficiencyArrayPixels, bins, maxMom)

	print "\n[PLOTTING GRAPHS]"
	graphs_2d("strips", efficiencyArrayStrips, bins, maxMom, args.strips)
	graphs_2d("pixels", efficiencyArrayPixels, bins, maxMom, args.pixels)
	plot_1D_graph("mom", [BCGraphMomStrips, ACGraphMomStrips, BCGraphMomPixels, ACGraphMomPixels], args.outputDirectory)
	plot_1D_graph("theta", [BCGraphThetaStrips, ACGraphThetaStrips, BCGraphThetaPixels, ACGraphThetaPixels],
				args.outputDirectory)
	res_graph(args.strips, momSigmaArrayStrips, momSigmaArrayPixels, bins, maxMom)

	print "[DONE]"

if __name__=='__main__':
	main()

