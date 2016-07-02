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
from collections import defaultdict

def get_MCParticle_details(MCparticle):
	MCMomTot = math.sqrt(math.pow(MCparticle.getMomentum()[0],2) + math.pow(MCparticle.getMomentum()[1],2) + math.pow(MCparticle.getMomentum()[2],2))
	#MCMomTrans = math.sqrt(math.pow(MCparticle.getMomentum()[0],2) + math.pow(MCparticle.getMomentum()[1],2))
	MCtheta = 90 - abs(math.degrees(math.atan((MCparticle.getMomentum()[2])/(math.sqrt(math.pow(MCparticle.getMomentum()[0],2) + math.pow(MCparticle.getMomentum()[1],2))))))
	return MCMomTot, MCtheta

def get_track_details(track): # CAN INCLUDE THIS AT LATER DATE TO LOOK AT QUALITY OF TRACKS AS BUNCHES INCREASE!!!! SHOULD DO THIS!!!!
	trackTransMom = math.fabs((((2.99792*math.pow(10,11))*(math.pow(10,-15))*5)/track.getOmega()))
	trackMom = trackTransMom*math.sqrt(1+pow(track.getTanLambda(),2))
	return trackTransMom

def tracking_analysis(MCparticleCollection, linkCollection, thetaCut, momCut):
	multiWeightsAbove = defaultdict(list)
	multiWeightsBelow = defaultdict(list)
	partialWeightsAbove = []
	partialWeightsBelow = []
	lessValues = [0,0,0,0,0] # totPions, totTracks, true, partrue, fake
	abvValues = [0,0,0,0,0]
	for MCparticle in MCparticleCollection:
		if ((MCparticle.getPDG() == 211 or MCparticle.getPDG() == -211)) and MCparticle.getGeneratorStatus() != 0:
			MCMomTot, MCtheta = get_MCParticle_details(MCparticle)
			if MCtheta > thetaCut and MCtheta < 86:
				linkCount = 0
				linkWeights = []
				for link in linkCollection:
					if link.getTo() == MCparticle:
						linkCount += 1
						linkWeights.append(link.getWeight())
				if linkCount == 1 and linkWeights[0] == 1: #If we find a pion and it has 1 track associated with it with weight of 1.
					if MCMomTot < momCut:
						lessValues[0] += 1
						lessValues[1] += 1
						lessValues[2] += 1
					if MCMomTot >= momCut:
						abvValues[0] += 1
						abvValues[1] += 1
						abvValues[2] += 1
				if linkCount == 1 and linkWeights[0] != 1: #If we find a pion and it has 1 track associated with it with weight not of 1.
					partialWeightsBelow.extend(linkWeights)
					if MCMomTot < momCut:
						lessValues[0] += 1
						lessValues[1] += 1
						lessValues[3] += 1
					if MCMomTot >= momCut:
						abvValues[0] += 1
						abvValues[1] += 1
						abvValues[3] += 1	
				if linkCount > 1: #If we find a pion and it has more than 1 tracks associated with it. These tracks are counted as fake.
					if MCMomTot < momCut:
						multiWeightsBelow["less_"+str(linkCount)].extend(linkWeights)
						lessValues[0] += 1
						lessValues[1] += linkCount
						lessValues[4] += linkCount
					if MCMomTot >= momCut:
						multiWeightsBelow["abov_"+str(linkCount)].extend(linkWeights)
						abvValues[0] += 1
						abvValues[1] += linkCount
						abvValues[4] += linkCount
				if linkCount == 0: #If we find a pion and it has no tracks associated with it. 
					if MCMomTot < momCut:
						lessValues[0] += 1
					if MCMomTot >= momCut:
						abvValues[0] += 1
						
			elif MCtheta >= 86:
				linkCount = 0
				linkWeights = []
				for link in linkCollection:
					if link.getTo() == MCparticle:
						linkCount += 1
						linkWeights.append(link.getWeight())
				if linkCount == 1 and linkWeights[0] != 1: #If we find a pion and it has 1 track associated with it with weight not of 1.
					partialWeightsAbove.extend(linkWeights)
				if linkCount > 1: #If we find a pion and it has more than 1 tracks associated with it. These tracks are counted as fake.
					if MCMomTot < momCut:
						multiWeightsAbove["less_"+str(linkCount)].extend(linkWeights)
					if MCMomTot >= momCut:
						multiWeightsAbove["abov_"+str(linkCount)].extend(linkWeights)


	return abvValues, lessValues, multiWeightsAbove, partialWeightsAbove, multiWeightsBelow, partialWeightsBelow

def track_validation(fileName): # COULD REWRITE THIS SO IT DOES NOT LOOP THROUGH EVERY PARTICLE, using "if in" i think!!!
	reader = IOIMPL.LCFactory.getInstance().createLCReader()
	reader.open(fileName)	
	values = [[0,0,0,0,0] for x in range(4)]# [[totPions, totTracks, true, partrue, fake] for abvEvent, lessEvent, abvOverlay, lessOverlay]
	multiWeightsAbove = defaultdict(list) 
	partialWeightsAbove = []
	multiWeightsBelow = defaultdict(list) 
	partialWeightsBelow = []
	for event in reader:
		MCParticlesSignal = event.getCollection("MCParticlesSkimmedEvent")
		MCParticlesOverlay = event.getCollection("MCParticlesSkimmedOverlay")
		truthLinkSignal = event.getCollection("TrackMCTruthLinkEvent")
		truthLinkOverlay = event.getCollection("TrackMCTruthLinkOverlay")
		# Checks to make sure things are like I think they should be!!!!
		trackNum = event.getCollection("Tracks").getNumberOfElements()
		if truthLinkSignal.getNumberOfElements() != trackNum or truthLinkOverlay.getNumberOfElements() != trackNum:
			print "Error: Different number of track elements to truthLink elements!!!"
			sys.exit(1)
		# Look at the event and overlay stuff seperately.
		abvValuesSignal, lessValuesSignal, multiWeightsSignalAbove, partialWeightsSignalAbove, multiWeightsSignalBelow, partialWeightsSignalBelow = tracking_analysis(MCParticlesSignal,truthLinkSignal,50,1)
		abvValuesOverlay, lessValuesOverlay, multiWeightsOverlayAbove, partialWeightsOverlayAbove, multiWeightsOverlayBelow, partialWeightsOverlayBelow = tracking_analysis(MCParticlesOverlay,truthLinkOverlay,50,1)
		partialWeightsAbove.extend(partialWeightsSignalAbove)
		partialWeightsAbove.extend(partialWeightsOverlayAbove)
		partialWeightsBelow.extend(partialWeightsSignalBelow)
		partialWeightsBelow.extend(partialWeightsOverlayBelow)
		for key in multiWeightsSignalAbove:
			multiWeightsAbove[key].extend(multiWeightsSignalAbove[key])
		for key in multiWeightsOverlayAbove:
			multiWeightsAbove[key].extend(multiWeightsOverlayAbove[key])
		for key in multiWeightsSignalBelow:
			multiWeightsBelow[key].extend(multiWeightsSignalBelow[key])
		for key in multiWeightsOverlayBelow:
			multiWeightsBelow[key].extend(multiWeightsOverlayBelow[key])
		temp = [abvValuesSignal,lessValuesSignal,abvValuesOverlay,lessValuesOverlay]
		for counter1 in range(4):
			for counter2 in range(5):
				values[counter1][counter2] += temp[counter1][counter2]
	return values, multiWeightsAbove, partialWeightsAbove, multiWeightsBelow, partialWeightsBelow

def parse_args():
	# Takes in arguments from the command line when script is called and returns them to main().
	currentDir = os.getcwd()
	parser = argparse.ArgumentParser(description='Processes .slcio file...'
						  ,epilog='In case of questions or problems, contact jt12194@my.bristol.ac.uk')

	parser.add_argument('-1', '--bunch1', help='Path to 1Bunch directory')

	parser.add_argument('-5', '--bunch5', help='Path to 5Bunch directory')

	parser.add_argument('-10', '--bunch10', help='Path to 10Bunch directory')

	parser.add_argument('-20', '--bunch20', help='Path to 20Bunch directory')

	parser.add_argument('-50', '--bunch50', help='Path to 50Bunch directory')

	return parser.parse_args()

def draw_graphs(data, bunchtypes):
	print "Making Graphs!!!"
	graphs = [ROOT.TGraphErrors() for x in range(18)]
	'''
	for index1 in range(6):
		for index2 in range(3):
			for index3 in range(2):
				for index4 in range(5):
	'''
	
	for counter in range(5):
		graphs[0].SetTitle(">1GeV, Efficiency")
		graphs[0].SetName(">1GeV, Efficiency")
		graphs[0].SetMinimum(0.85)
		graphs[0].SetMaximum(1)
		graphs[0].SetPoint(counter, bunchtypes[counter], data[counter][0][0][0])
		graphs[0].SetPointError(counter, 0, data[counter][0][0][1]) 
		graphs[1].SetTitle(">1GeV, Partial Rate")
		graphs[1].SetName(">1GeV, Partial Rate")
		graphs[1].SetMinimum(pow(10,-3))
		graphs[1].SetMaximum(1)
		graphs[1].SetPoint(counter, bunchtypes[counter], data[counter][0][0][2])
		graphs[1].SetPointError(counter, 0, data[counter][0][0][3])
		graphs[2].SetTitle(">1GeV, MultiRate")
		graphs[2].SetName(">1GeV, MultiRate")
		graphs[2].SetMinimum(pow(10,-5))
		graphs[2].SetMaximum(pow(10,-2))
		graphs[2].SetPoint(counter, bunchtypes[counter], data[counter][0][0][4])
		graphs[2].SetPointError(counter, 0, data[counter][0][0][5])
		graphs[3].SetTitle("Above1GeVEfficiencyOverlay")
		graphs[3].SetName("Above1GeVEfficiencyOverlay")
		graphs[3].SetMinimum(0)
		graphs[3].SetMaximum(1)
		graphs[3].SetPoint(counter, bunchtypes[counter], data[counter][1][0][0])
		graphs[3].SetPointError(counter, 0, data[counter][1][0][1])
		graphs[4].SetTitle("Above1GeVPartialOverlay")
		graphs[4].SetName("Above1GeVPartialOverlay")
		graphs[4].SetMinimum(pow(10,-3))
		graphs[4].SetMaximum(1)
		graphs[4].SetPoint(counter, bunchtypes[counter], data[counter][1][0][2])
		graphs[4].SetPointError(counter, 0, data[counter][1][0][3])
		graphs[5].SetTitle("Above1GeVFakeOverlay")
		graphs[5].SetName("Above1GeVFakeOverlay")
		graphs[5].SetMinimum(pow(10,-5))
		graphs[5].SetMaximum(pow(10,-2))
		graphs[5].SetPoint(counter, bunchtypes[counter], data[counter][1][0][4])
		graphs[5].SetPointError(counter, 0, data[counter][1][0][5])
		graphs[6].SetTitle("Above1GeVEfficiencyTot")
		graphs[6].SetName("Above1GeVEfficiencyTot")
		graphs[6].SetMinimum(0)
		graphs[6].SetMaximum(1)
		graphs[6].SetPoint(counter, bunchtypes[counter], data[counter][2][0][0])
		graphs[6].SetPointError(counter, 0, data[counter][2][0][1])
		graphs[7].SetTitle("Above1GeVPartialTot")
		graphs[7].SetName("Above1GeVPartialTot")
		graphs[7].SetMinimum(pow(10,-3))
		graphs[7].SetMaximum(1)
		graphs[7].SetPoint(counter, bunchtypes[counter], data[counter][2][0][2])
		graphs[7].SetPointError(counter, 0, data[counter][2][0][3])
		graphs[8].SetTitle("Above1GeVFakeTot")
		graphs[8].SetName("Above1GeVFakeTot")
		graphs[8].SetMinimum(pow(10,-5))
		graphs[8].SetMaximum(pow(10,-2))
		graphs[8].SetPoint(counter, bunchtypes[counter], data[counter][2][0][4])
		graphs[8].SetPointError(counter, 0, data[counter][2][0][5])
		graphs[9].SetTitle("<1GeV, Efficiency")
		graphs[9].SetName("<1GeV, Efficiency")
		graphs[9].SetMinimum(0.4)
		graphs[9].SetMaximum(0.8)
		graphs[9].SetPoint(counter, bunchtypes[counter], data[counter][0][1][0])
		graphs[9].SetPointError(counter, 0, data[counter][0][1][1])
		graphs[10].SetTitle("<1GeV, Partial Rate")
		graphs[10].SetName("<1GeV, Partial Rate")
		graphs[10].SetMinimum(pow(10,-3))
		graphs[10].SetMaximum(1)
		graphs[10].SetPoint(counter, bunchtypes[counter], data[counter][0][1][2])
		graphs[10].SetPointError(counter, 0, data[counter][0][1][3])
		graphs[11].SetTitle("<1GeV, MultiRate")
		graphs[11].SetName("<1GeV, MultiRate")
		graphs[11].SetMinimum(pow(10,-5))
		graphs[11].SetMaximum(pow(10,-2))
		graphs[11].SetPoint(counter, bunchtypes[counter], data[counter][0][1][4])
		graphs[11].SetPointError(counter, 0, data[counter][0][1][5])
		graphs[12].SetTitle("Less1GeVEfficiencyOverlay")
		graphs[12].SetName("Less1GeVEfficiencyOverlay")
		graphs[12].SetMinimum(0)
		graphs[12].SetMaximum(1)
		graphs[12].SetPoint(counter, bunchtypes[counter], data[counter][1][1][0])
		graphs[12].SetPointError(counter, 0, data[counter][1][1][1])
		graphs[13].SetTitle("Less1GeVPartialOverlay")
		graphs[13].SetName("Less1GeVPartialOverlay")
		graphs[13].SetMinimum(pow(10,-3))
		graphs[13].SetMaximum(1)
		graphs[13].SetPoint(counter, bunchtypes[counter], data[counter][1][1][2])
		graphs[13].SetPointError(counter, 0, data[counter][1][1][3])
		graphs[14].SetTitle("Less1GeVFakeOverlay")
		graphs[14].SetName("Less1GeVFakeOverlay")
		graphs[14].SetMinimum(pow(10,-5))
		graphs[14].SetMaximum(pow(10,-2))
		graphs[14].SetPoint(counter, bunchtypes[counter], data[counter][1][1][4])
		graphs[14].SetPointError(counter, 0, data[counter][1][1][5])
		graphs[15].SetTitle("Less1GeVEfficiencyTot")
		graphs[15].SetName("Less1GeVEfficiencyTot")
		graphs[15].SetMinimum(0)
		graphs[15].SetMaximum(1)
		graphs[15].SetPoint(counter, bunchtypes[counter], data[counter][2][1][0])
		graphs[15].SetPointError(counter, 0, data[counter][2][1][1])
		graphs[16].SetTitle("Less1GeVPartialTot")
		graphs[16].SetName("Less1GeVPartialTot")
		graphs[16].SetMinimum(pow(10,-3))
		graphs[16].SetMaximum(1)
		graphs[16].SetPoint(counter, bunchtypes[counter], data[counter][2][1][2])
		graphs[16].SetPointError(counter, 0, data[counter][2][1][3])
		graphs[17].SetTitle("Less1GeVFakeTot")
		graphs[17].SetName("Less1GeVFakeTot")
		graphs[17].SetMinimum(pow(10,-5))
		graphs[17].SetMaximum(pow(10,-2))
		graphs[17].SetPoint(counter, bunchtypes[counter], data[counter][2][1][4])
		graphs[17].SetPointError(counter, 0, data[counter][2][1][5])

	for graph in graphs:
		#graph.GetXaxis().SetTitle("Number of Bunches")
		graph.GetYaxis().SetTitle("Rate")
		graph.GetXaxis().SetTitle("Number of Bunches")
		graph.GetXaxis().SetLimits(0, 51)
	'''
	sets = [[0,3,6],[9,12,15],[1,4,7],[10,13,16],[2,5,8],[11,14,17]]
	graphNames = ["abvEfficiency", "lessEfficiency", "abvPartial", "lessPartial", "abvFake", "lessFake"]
	
	for count1 in range(3):
		c = TCanvas("C", "Canvas", 500, 1000)
		c.Divide(1,2)
		c.cd(1)
		for count2 in range(3):
			if count2 == 0:
				graphs[sets[count1*2][count2]].SetMarkerColor(count2+2)
				graphs[sets[count1*2][count2]].SetLineColor(count2+2)
				graphs[sets[count1*2][count2]].SetFillColor(count2+2)
				graphs[sets[count1*2][count2]].Draw("ALP")
			else:
				graphs[sets[count1*2][count2]].SetMarkerColor(count2+2)
				graphs[sets[count1*2][count2]].SetLineColor(count2+2)
				graphs[sets[count1*2][count2]].SetFillColor(count2+2)
				graphs[sets[count1*2][count2]].Draw("SAME")
		c.cd(2)
		for count2 in range(3):
			if count2 == 0:
				graphs[sets[(count1*2)+1][count2]].SetMarkerColor(count2+2)
				graphs[sets[(count1*2)+1][count2]].SetLineColor(count2+2)
				graphs[sets[(count1*2)+1][count2]].SetFillColor(count2+2)
				graphs[sets[(count1*2)+1][count2]].Draw("ALP")
			else:
				graphs[sets[(count1*2)+1][count2]].SetMarkerColor(count2+2)
				graphs[sets[(count1*2)+1][count2]].SetLineColor(count2+2)
				graphs[sets[(count1*2)+1][count2]].SetFillColor(count2+2)
				graphs[sets[(count1*2)+1][count2]].Draw("SAME")
		
		c.Modified()
		c.Update()
		raw_input()
		img = ROOT.TImage.Create()
		img.FromPad(c)
		img.WriteImage(graphNames[count1]+".png")
		del c, img
	'''
	fits(graphs[0])
	fits(graphs[9])
	fits(graphs[1])
	fits(graphs[10])
	fits(graphs[2])
	fits(graphs[11])
	

def fits(graph):
	print "For graph name -> " + graph.GetName()
	pol1Fit = ROOT.TF1("pol1Fit", "pol1", 0, 50)
	expoFit = ROOT.TF1("expoFit", "expo", 0, 50)
	graph.Fit(pol1Fit, "RNQF")	
	print "pol1Chi -> " + str(pol1Fit.GetChisquare())
	graph.Fit(expoFit, "RNQF")
	print "expoChi -> " + str(expoFit.GetChisquare())
	# Now need to create the extrapolation graphs with errors out to 2450!!!
	extPol1 = ROOT.TGraphErrors()
	extPol1.SetName("extpol1")
	extExpo = TGraphAsymmErrors()
	extExpo.SetName("extExpo")
	for point in range(20):
		xCoord = ((2450/19)*point)
		#pol1
		extPol1.SetPoint(point, xCoord, pol1Fit.GetParameter(0) + (pol1Fit.GetParameter(1)*xCoord))
		extPol1.SetPointError(point, 0, pol1Fit.GetParError(0) + (pol1Fit.GetParError(1)*xCoord))
		#expo
		lowError = expoFit.Eval(xCoord) - math.exp((expoFit.GetParameter(0)-expoFit.GetParError(0))+((expoFit.GetParameter(1)-expoFit.GetParError(1))*xCoord))
		upError = math.exp((expoFit.GetParameter(0)+expoFit.GetParError(0))+((expoFit.GetParameter(1)+expoFit.GetParError(1))*xCoord)) - expoFit.Eval(xCoord)
		extExpo.SetPoint(point, xCoord, math.exp(expoFit.GetParameter(0) + (expoFit.GetParameter(1)*xCoord)))
		extExpo.SetPointError(point, 0, 0, lowError, upError)

	c1 = TCanvas("C", "Canvas")
	c1.SetGrid()	
	c1.SetTitle(graph.GetName())
	extPol1.GetYaxis().SetTitle("Rate")
	extPol1.GetXaxis().SetTitle("Number of Bunches")
	extPol1.SetFillColor(2)
	extPol1.SetFillStyle(3006)
	extExpo.SetFillColor(4)
	extExpo.SetFillStyle(3003)
	extPol1.SetTitle(graph.GetName())
	extExpo.SetTitle(graph.GetName())

	extExpo.Draw("aL3")
	extPol1.Draw("SAME L3")
	graph.Draw("SAME L")

	c1.Modified()
	c1.Update()
	raw_input()
	img = ROOT.TImage.Create()
	img.FromPad(c1)
	img.WriteImage(graph.GetName()+".png")
	del c1, pol1Fit, expoFit, extPol1, extExpo


def draw_weights(name, weights, bins, minX, maxX):
	theHistAbv = ROOT.TH2F(name+"Abv", name+"Abv", 50, 0, 1, 10, 0, 10)
	theHistLess = ROOT.TH2F(name+"Less", name+"Less", 50, 0, 1, 10, 0, 10)
	for key in weights:
		size = len(weights[key])
		num = int(key[-1:])
		energyType = key[:4]
		for dataPoint in weights[key]:
			if energyType == "abov":
				theHistAbv.Fill(dataPoint, num-0.5)	
			if energyType == "less":
				theHistLess.Fill(dataPoint, num-0.5)

	c = TCanvas("c", name)
	theHistAbv.Draw("LEGO")
	theHistAbv.GetXaxis().SetTitle("weight")
	theHistAbv.GetYaxis().SetTitle("Number of Multitracks")
	raw_input()
	img = ROOT.TImage.Create()
	img.FromPad(c)
	img.WriteImage(name +"Abv.png")
	del c

	c = TCanvas("c", name)
	theHistLess.Draw("LEGO")
	theHistLess.GetXaxis().SetTitle("weight")
	theHistLess.GetYaxis().SetTitle("Number of Multitracks")
	raw_input()
	img = ROOT.TImage.Create()
	img.FromPad(c)
	img.WriteImage(name +"Less.png")
	del c

def plot_simple_histogram(name, data, bins, minX, maxX, fit, centreFit, rangeFit):
	print "Plotting-> " + name 
	# Plots histogram and fits anything according to arguments.
	hist = ROOT.TH1F(name, name, bins, minX, maxX)
	c = TCanvas("c", name)
	for dataPoint in data:
		hist.Fill(dataPoint)

	hist.GetXaxis().SetTitle("weight")
	hist.GetYaxis().SetTitle("Frequency")
	
	hist.Draw()
	hist.Print()

	c.Modified()
	c.Update()
	c.SetGrid()

	if fit:
		hist.Fit("gaus","","",centreFit - rangeFit, centreFit + rangeFit)

	c.Modified()
	c.Update()
	hist.SetStats(False)

	raw_input()
	img = ROOT.TImage.Create()
	img.FromPad(c)
	img.WriteImage(name + ".png")

def main():
	args = parse_args()
	if not args:
		print 'Invalid Arguments'
		sys.exit(1)

	files = [[] for x in range(5)]
	files[0] = commonFunctions.input_files(args.bunch1, ".slcio")
	files[1] = commonFunctions.input_files(args.bunch5, ".slcio")
	files[2] = commonFunctions.input_files(args.bunch10, ".slcio")
	files[3] = commonFunctions.input_files(args.bunch20, ".slcio")
	files[4] = commonFunctions.input_files(args.bunch50, ".slcio")

	bunchtypes = [1, 5, 10, 20, 50]
	count = 0
	data = [[[[0,0,0,0,0,0] for x in range(2)] for x in range(3)] for x in range(5)] #[[[abv,less]signal,overlay,tot]bunch1,bunch5 etc....]
	multiWeightsAbove = defaultdict(list)
	partialWeightsAbove = []
	multiWeightsBelow = defaultdict(list)
	partialWeightsBelow = []
	temp = [0 for x in range(18)]
	for bunchType in files:
		print "Analysing Bunch type -> " + str(bunchtypes[count])
		values = [[0,0,0,0,0] for x in range(4)] 
		fileCount = 1
		for fileName in bunchType:
			print "Processing File num -> " + str(fileCount)
			fileCount += 1
			valuesTemp, multiWeightsAboveTemp, partialWeightsAboveTemp, multiWeightsBelowTemp, partialWeightsBelowTemp = track_validation(fileName)
			partialWeightsAbove.extend(partialWeightsAboveTemp)
			for key in multiWeightsAboveTemp:
				multiWeightsAbove[key].extend(multiWeightsAboveTemp[key])
			partialWeightsBelow.extend(partialWeightsBelowTemp)
			for key in multiWeightsBelowTemp:
				multiWeightsBelow[key].extend(multiWeightsBelowTemp[key])
			for counter1 in range(4):
				for counter2 in range(5):
					values[counter1][counter2] += valuesTemp[counter1][counter2]

		print "Calculating data Values..."
		# values = [[0,0,0,0,0] for x in range(4)]
		# data[bunches][type][cut][value] # values[type][value]
		# Above,signal
		print "Total pions -> " +str(values[0][2])
		data[count][0][0][0] = values[0][2] / values[0][0] # Efficiency
		data[count][0][0][2] = values[0][3] / values[0][0] # Efficiency		
		data[count][0][0][4] = values[0][4] / values[0][1] # Efficiency
		
		# less,signal
		data[count][0][1][0] = values[1][2] / values[1][0] # Efficiency		
		data[count][0][1][2] = values[1][3] / values[1][0] # Efficiency		
		data[count][0][1][4] = values[1][4] / values[1][1] # Efficiency]
		
		# Above,overlay
		data[count][1][0][0] = values[2][2] / values[2][0] # Efficiency	
		data[count][1][0][2] = values[2][3] / values[2][0] # Efficiency	
		data[count][1][0][4] = values[2][4] / values[2][1] # Efficiency
		
		# less,overlay
		data[count][1][1][0] = values[3][2] / values[3][0] # Efficiency		
		data[count][1][1][2] = values[3][3] / values[3][0] # Efficiency	
		data[count][1][1][4] = values[3][4] / values[3][1] # Efficiency
		
		# Above,tot
		data[count][2][0][0] = (values[0][2]+values[2][2]) / (values[0][0]+values[2][0]) # Efficiency	
		data[count][2][0][2] = (values[0][3]+values[2][3]) / (values[0][0]+values[2][0]) # Efficiency	
		data[count][2][0][4] = (values[0][4]+values[2][4]) / (values[0][1]+values[2][1]) # Efficiency
		
		# less,tot
		data[count][2][1][0] = (values[1][2]+values[3][2]) / (values[1][0]+values[3][0]) # Efficiency		
		data[count][2][1][2] = (values[1][3]+values[3][3]) / (values[1][0]+values[3][0]) # Efficiency	
		data[count][2][1][4] = (values[1][4]+values[3][4]) / (values[1][1]+values[3][1]) # Efficiency

		#OLD ONES
		
		'''
		data[count][0][0][1] = math.sqrt(((values[0][2] / values[0][0])*(1-(values[0][2] / values[0][0])))/(values[0][0]))
		data[count][0][0][3] = math.sqrt(((values[0][3] / values[0][0])*(1-(values[0][3] / values[0][0])))/(values[0][0]))
		data[count][0][0][5] = math.sqrt(((values[0][4] / values[0][1])*(1-(values[0][4] / values[0][1])))/(values[0][1]))
		data[count][0][1][1] = math.sqrt(((values[1][2] / values[1][0])*(1-(values[1][2] / values[1][0])))/(values[1][0]))
		data[count][0][1][3] = math.sqrt(((values[1][3] / values[1][0])*(1-(values[1][3] / values[1][0])))/(values[1][0]))
		data[count][0][1][5] = math.sqrt(((values[1][4] / values[1][1])*(1-(values[1][4] / values[1][1])))/(values[1][1]))
		data[count][1][0][1] = math.sqrt(((values[2][2] / values[2][0])*(1-(values[2][2] / values[2][0])))/(values[2][0]))
		data[count][1][0][3] = math.sqrt(((values[2][3] / values[2][0])*(1-(values[2][3] / values[2][0])))/(values[2][0]))
		data[count][1][0][5] = math.sqrt(((values[2][4] / values[2][1])*(1-(values[2][4] / values[2][1])))/(values[2][1]))
		data[count][1][1][1] = math.sqrt(((values[3][2] / values[3][0])*(1-(values[3][2] / values[3][0])))/(values[3][0]))
		data[count][1][1][3] = math.sqrt(((values[3][3] / values[3][0])*(1-(values[3][3] / values[3][0])))/(values[3][0]))
		data[count][1][1][5] = math.sqrt(((values[3][4] / values[3][1])*(1-(values[3][4] / values[3][1])))/(values[3][1]))
		data[count][2][0][1] = math.sqrt((((values[0][2]+values[2][2]) / (values[0][0]+values[2][0]))*(1-((values[0][2]+values[2][2]) / (values[0][0]+values[2][0]))))/((values[0][0]+values[2][0])))
		data[count][2][0][3] = math.sqrt((((values[0][3]+values[2][3]) / (values[0][0]+values[2][0]))*(1-((values[0][3]+values[2][3]) / (values[0][0]+values[2][0]))))/((values[0][0]+values[2][0])))
		data[count][2][0][5] = math.sqrt((((values[0][4]+values[2][4]) / (values[0][1]+values[2][1]))*(1-((values[0][4]+values[2][4]) / (values[0][1]+values[2][1]))))/((values[0][1]+values[2][1])))
		data[count][2][1][1] = math.sqrt((((values[1][2]+values[3][2]) / (values[1][0]+values[3][0]))*(1-((values[1][2]+values[3][2]) / (values[1][0]+values[3][0]))))/((values[1][0]+values[3][0])))
		data[count][2][1][3] = math.sqrt((((values[1][3]+values[3][3]) / (values[1][0]+values[3][0]))*(1-((values[1][3]+values[3][3]) / (values[1][0]+values[3][0]))))/((values[1][0]+values[3][0])))
		data[count][2][1][5] = math.sqrt((((values[1][4]+values[3][4]) / (values[1][1]+values[3][1]))*(1-((values[1][4]+values[3][4]) / (values[1][1]+values[3][1]))))/((values[1][1]+values[3][1])))
		'''
		#NEWEST ONES!!!!!
		
		data[count][0][0][1] = math.sqrt(((math.fabs(values[0][2]-temp[0]) / values[0][0])*(1-(math.fabs(values[0][2]-temp[0]) / values[0][0])))/(values[0][0]))
		data[count][0][0][3] = math.sqrt(((math.fabs(values[0][3]-temp[1]) / values[0][0])*(1-(math.fabs(values[0][3]-temp[1]) / values[0][0])))/(values[0][0]))
		data[count][0][0][5] = math.sqrt(((math.fabs(values[0][4]-temp[2]) / values[0][1])*(1-(math.fabs(values[0][4]-temp[2]) / values[0][1])))/(values[0][1]))
		data[count][0][1][1] = math.sqrt(((math.fabs(values[1][2]-temp[3]) / values[1][0])*(1-(math.fabs(values[1][2]-temp[3]) / values[1][0])))/(values[1][0]))
		data[count][0][1][3] = math.sqrt(((math.fabs(values[1][3]-temp[4]) / values[1][0])*(1-(math.fabs(values[1][3]-temp[4]) / values[1][0])))/(values[1][0]))
		data[count][0][1][5] = math.sqrt(((math.fabs(values[1][4]-temp[5]) / values[1][1])*(1-(math.fabs(values[1][4]-temp[5]) / values[1][1])))/(values[1][1]))
		data[count][1][0][1] = math.sqrt(((math.fabs(values[2][2]-temp[6]) / values[2][0])*(1-(math.fabs(values[2][2]-temp[6]) / values[2][0])))/(values[2][0]))
		data[count][1][0][3] = math.sqrt(((math.fabs(values[2][3]-temp[7]) / values[2][0])*(1-(math.fabs(values[2][3]-temp[7]) / values[2][0])))/(values[2][0]))
		data[count][1][0][5] = math.sqrt(((math.fabs(values[2][4]-temp[8]) / values[2][1])*(1-(math.fabs(values[2][4]-temp[8]) / values[2][1])))/(values[2][1]))
		data[count][1][1][1] = math.sqrt(((math.fabs(values[3][2]-temp[9]) / values[3][0])*(1-(math.fabs(values[3][2]-temp[9])/ values[3][0])))/(values[3][0]))
		data[count][1][1][3] = math.sqrt(((math.fabs(values[3][3]-temp[10]) / values[3][0])*(1-(math.fabs(values[3][3]-temp[10]) / values[3][0])))/(values[3][0]))
		data[count][1][1][5] = math.sqrt(((math.fabs(values[3][4]-temp[11]) / values[3][1])*(1-(math.fabs(values[3][4]-temp[11]) / values[3][1])))/(values[3][1]))
		data[count][2][0][1] = math.sqrt(((math.fabs((values[0][2]+values[2][2])-temp[12]) / (values[0][0]+values[2][0]))*(1-(math.fabs((values[0][2]+values[2][2])-temp[12]) / (values[0][0]+values[2][0]))))/((values[0][0]+values[2][0])))
		data[count][2][0][3] = math.sqrt(((math.fabs((values[0][3]+values[2][3])-temp[13]) / (values[0][0]+values[2][0]))*(1-(math.fabs((values[0][3]+values[2][3])-temp[13]) / (values[0][0]+values[2][0]))))/((values[0][0]+values[2][0])))
		data[count][2][0][5] = math.sqrt(((math.fabs((values[0][4]+values[2][4])-temp[14]) / (values[0][1]+values[2][1]))*(1-(math.fabs((values[0][4]+values[2][4])-temp[14]) / (values[0][1]+values[2][1]))))/((values[0][1]+values[2][1])))
		data[count][2][1][1] = math.sqrt(((math.fabs((values[1][2]+values[3][2])-temp[15]) / (values[1][0]+values[3][0]))*(1-(math.fabs((values[1][2]+values[3][2])-temp[15]) / (values[1][0]+values[3][0]))))/((values[1][0]+values[3][0])))
		data[count][2][1][3] = math.sqrt(((math.fabs((values[1][3]+values[3][3])-temp[16]) / (values[1][0]+values[3][0]))*(1-(math.fabs((values[1][3]+values[3][3])-temp[16]) / (values[1][0]+values[3][0]))))/((values[1][0]+values[3][0])))
		data[count][2][1][5] = math.sqrt(((math.fabs((values[1][4]+values[3][4])-temp[17]) / (values[1][1]+values[3][1]))*(1-(math.fabs((values[1][4]+values[3][4])-temp[17]) / (values[1][1]+values[3][1]))))/((values[1][1]+values[3][1])))
		

		temp[0] = values[0][2]
		temp[1] = values[0][3]
		temp[2] = values[0][4]
		temp[3] = values[1][2]
		temp[4] = values[1][3]
		temp[5] = values[1][4]
		temp[6] = values[2][2]
		temp[7] = values[2][3]
		temp[8] = values[2][4]
		temp[9] = values[3][2]
		temp[10] = values[3][3]
		temp[11] = values[3][4]
		temp[12] = values[0][2]+values[2][2]
		temp[13] = values[0][3]+values[2][3]
		temp[14] = values[0][4]+values[2][4]
		temp[15] = values[1][2]+values[3][2]
		temp[16] = values[1][3]+values[3][3]
		temp[17] = values[1][4]+values[3][4]
		
		count += 1


	#plot_simple_histogram("partialWeightsAbove", partialWeightsAbove, 100, 0.5, 1, False, 0, 0)
	plot_simple_histogram("partialWeightsBelow", partialWeightsBelow, 100, 0.5, 1, False, 0, 0)
	#draw_weights("multiWeightsAbove", multiWeightsAbove, 100, 0, 1)
	#draw_weights("multiWeightsBelow", multiWeightsBelow, 100, 0, 1)
	draw_graphs(data, bunchtypes)

if __name__=='__main__':
	main()