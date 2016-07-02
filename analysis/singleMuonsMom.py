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
from collections import defaultdict

def get_data(fileName): #COMBINE WITH ABOVE FUNCTION!!!!!!
	reader = IOIMPL.LCFactory.getInstance().createLCReader()
	reader.open(fileName)

	eventCounter = 0
	deltaMom = []

	for event in reader:
		MCParticles = event.getCollection("MCParticle")
		track = event.getCollection("Tracks")
		muonCount = 0
		trackCount = 0
		for particle in MCParticles:
			if particle.getPDG() == 13 or particle.getPDG() == -13:
				muonCount += 1
				mom = particle.getMomentum()
				totalMom = math.sqrt(math.pow(mom[0],2) + math.pow(mom[1],2) + math.pow(mom[2],2))
				transMom = math.sqrt(math.pow(mom[0],2) + math.pow(mom[1],2))
				theta = round(math.degrees(math.acos(mom[2]/totalMom))) 
				if eventCounter == 0:
					theta = round(math.degrees(math.acos(mom[2]/totalMom)))
					print "Total Momentum = " + str(totalMom) + "\nTransverse Momentum = " + str(transMom) + "\nTheta = " + str(theta)

		for particleTrack in track:
			trackCount += 1
			omega = particleTrack.getOmega()
			trackTransMom = math.fabs((((2.99792*math.pow(10,11))*(math.pow(10,-15))*5)/omega)) 

		if muonCount == 1 and trackCount == 1:
			deltaMom.append(trackTransMom - transMom)
			eventCounter += 1
		'''
		if muonCount == 0 or trackCount == 0:
			print "----Could not find a particle----"

		if muonCount > 1 or trackCount > 1:
			print "----More than 1 muon----"
		'''
		eventCounter += 1

	reader.close()

	print "Number of events used->" + str(eventCounter)

	return totalMom, transMom, theta, deltaMom

def comparison_graph(moms, values, directory):
	c = TCanvas()
	c.SetLogy()
	c.SetLogx()

	graphs = {}
	graphCount = 0
	leg = ROOT.TLegend()

	for description in moms:
		print description
		momArray = array.array("f", moms[description])
		valueArray = array.array("f",values[description])
		graphs[description] = ROOT.TGraph(7, momArray, valueArray)
		graphs[description].SetMarkerColor(graphCount + 1)
		graphs[description].SetName(description)
		leg.AddEntry(graphs[description],description,"p")
		if graphCount == 0:
			graphs[description].GetYaxis().SetTitle("#sigma(p_{T})/p_{t}^{2} [GeV^{-1}]")
			graphs[description].GetXaxis().SetTitle("p [GeV]")
			graphs[description].SetMinimum(pow(10,-5))
			graphs[description].SetMaximum(pow(10,-1))
			graphs[description].GetXaxis().SetLimits(0.9, 101)
			graphs[description].Draw("A*")
		else: graphs[description].Draw("*")

		c.Update()
		graphCount += 1

	leg.Draw()
	c.Update()

	imp = raw_input()
	img = ROOT.TImage.Create()
	img.FromPad(c)
	img.WriteImage(directory + "/comparison.png")

def mom_graph(delta_mom, outputname, bins, width):
	hist = ROOT.TH1F(outputname, outputname, bins, -width, width)
	for mom in delta_mom:
		hist.Fill(mom)
	hist.Fit("gaus")
	# print hist.GetFunction("gaus").GetParameter(2)
	# print hist.GetRMS()
	return hist.GetFunction("gaus").GetParameter(2)

def main():
	# main
	args = commonFunctions.simple_parse_args()

	if not args:
		print 'Invalid Arguments'
		sys.exit(1)

	fileList = commonFunctions.input_files(args.inputDir, ".slcio")

	values = defaultdict(list)
	moms = defaultdict(list)
	fileCount = 0

	for fileName in fileList:
		fileCount += 1
		print "File -> " + str(fileCount)

		outputName, eventCount, detector, energy = commonFunctions.get_fileData(fileName)

		name, extension = os.path.splitext(fileName)
		checkname = name[-25:]
		print checkname

		totalMom, transMom, theta, deltaMom = get_data(fileName)

		width = totalMom / 10
		value = mom_graph(deltaMom, outputName, args.numberOfBins, width)/math.pow(transMom,2)

		#moms[str(theta) + "_" + str(detector)].append(totalMom)
		#values[str(theta) + "_" + str(detector)].append(value)

		moms[str(theta) + "_" + checkname].append(totalMom)
		print str(theta) + "_" + checkname
		values[str(theta) + "_" + checkname].append(value)


	comparison_graph(moms, values, args.inputDir)

	print "[DONE]"

if __name__=='__main__':
	main()