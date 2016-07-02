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

def get_MCParticle_details(MCparticle):
	MCMomTot = math.sqrt(math.pow(MCparticle.getMomentum()[0],2) + math.pow(MCparticle.getMomentum()[1],2) + math.pow(MCparticle.getMomentum()[2],2))
	#MCMomTrans = math.sqrt(math.pow(MCparticle.getMomentum()[0],2) + math.pow(MCparticle.getMomentum()[1],2))
	MCtheta = 90 - abs(math.degrees(math.atan((MCparticle.getMomentum()[2])/(math.sqrt(math.pow(MCparticle.getMomentum()[0],2) + math.pow(MCparticle.getMomentum()[1],2))))))
	return MCMomTot, MCtheta

def get_track_details(track): # CAN INCLUDE THIS AT LATER DATE TO LOOK AT QUALITY OF TRACKS AS BUNCHES INCREASE!!!! SHOULD DO THIS!!!!
	trackTransMom = math.fabs((((2.99792*math.pow(10,11))*(math.pow(10,-15))*5)/track.getOmega()))
	trackMom = trackTransMom*math.sqrt(1+pow(track.getTanLambda(),2))
	return trackTransMom

def track_validation(fileName): # COULD REWRITE THIS SO IT DOES NOT LOOP THROUGH EVERY PARTICLE, using "if in" i think!!!
	reader = IOIMPL.LCFactory.getInstance().createLCReader()
	reader.open(fileName)	
	lessValues = [0,0,0,0,0] # totPions, totTracks, true, partrue, fake
	abvValues = [0,0,0,0,0]
	for event in reader:
		MCParticles = event.getCollection("MCParticlesSkimmed")
		truthLink = event.getCollection("TrackMCTruthLink")
		if truthLink.getNumberOfElements() != event.getCollection("Tracks").getNumberOfElements():
			print "Error: Different number of track elements to truthLink elements!!!"
			sys.exit(1)
		for MCparticle in MCParticles:
			if (MCparticle.getPDG() == 211 or MCparticle.getPDG() == -211):
				MCMomTot, MCtheta = get_MCParticle_details(MCparticle)
				if MCtheta > 50:
				#if MCtheta > 15 and MCtheta < 45:
					linkCount = 0
					for link in truthLink:
						if link.getTo() == MCparticle:
							linkCount += 1
							linkWeight = link.getWeight()

					if linkCount == 1 and linkWeight == 1: #If we find a pion and it has 1 track associated with it with weight of 1.
						if MCMomTot < 1:
							lessValues[0] += 1
							lessValues[1] += 1
							lessValues[2] += 1
						if MCMomTot >= 1:
							abvValues[0] += 1
							abvValues[1] += 1
							abvValues[2] += 1

					if linkCount == 1 and linkWeight != 1: #If we find a pion and it has 1 track associated with it with weight not of 1.
						if MCMomTot < 1:
							lessValues[0] += 1
							lessValues[1] += 1
							lessValues[3] += 1
						if MCMomTot >= 1:
							abvValues[0] += 1
							abvValues[1] += 1
							abvValues[3] += 1	

					if linkCount > 1: #If we find a pion and it has more than 1 tracks associated with it. These tracks are counted as fake.
						if MCMomTot < 1:
							lessValues[0] += 1
							lessValues[1] += linkCount
							lessValues[4] += linkCount
						if MCMomTot >= 1:
							abvValues[0] += 1
							abvValues[1] += linkCount
							abvValues[4] += linkCount

					if linkCount == 0: #If we find a pion and it has no tracks associated with it. 
						if MCMomTot < 1:
							lessValues[0] += 1
						if MCMomTot >= 1:
							abvValues[0] += 1


	return lessValues, abvValues

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

def draw_graphs(lessData, abvData, bunchtypes):
	print "Making Graphs!!!"
	graphs = [ROOT.TGraphErrors() for x in range(6)]

	for counter in range(5):
		graphs[0].SetName("<1GeVEfficiency")
		graphs[0].SetMinimum(0.4)
		graphs[0].SetMaximum(0.7)
		graphs[0].SetPoint(counter, bunchtypes[counter], lessData[counter][0]) # lessEff
		graphs[0].SetPointError(counter, 0, lessData[counter][1]) 
		graphs[1].SetName("<1GeVPartial")
		graphs[1].SetMinimum(0)
		graphs[1].SetMaximum(0.3)
		graphs[1].SetPoint(counter, bunchtypes[counter], lessData[counter][2]) # lessPar
		graphs[1].SetPointError(counter, 0, lessData[counter][3]) 
		graphs[2].SetName("<1GeVFake")
		graphs[2].SetMinimum(0)
		graphs[2].SetMaximum(0.1)
		graphs[2].SetPoint(counter, bunchtypes[counter], lessData[counter][4]) # lessFake
		graphs[2].SetPointError(counter, 0, lessData[counter][5]) 
		graphs[3].SetName(">=1GeVEfficiency")
		graphs[3].SetMinimum(0.8)
		graphs[3].SetMaximum(1)
		graphs[3].SetPoint(counter, bunchtypes[counter], abvData[counter][0]) # abvEff
		graphs[3].SetPointError(counter, 0, lessData[counter][1]) 
		graphs[4].SetName(">=1GeVPartial")
		graphs[4].SetMinimum(0)
		graphs[4].SetMaximum(0.1)
		graphs[4].SetPoint(counter, bunchtypes[counter], abvData[counter][2]) # abvPar
		graphs[4].SetPointError(counter, 0, lessData[counter][3]) 
		graphs[5].SetName(">=1GeVFake")
		graphs[5].SetMinimum(0)
		graphs[5].SetMaximum(0.1)
		graphs[5].SetPoint(counter, bunchtypes[counter], abvData[counter][4]) # abvFake
		graphs[5].SetPointError(counter, 0, lessData[counter][5]) 	

	fits = [ROOT.TF1("pol2Fit", "pol2", 0, 50),ROOT.TF1("expoFit", "expo", 0, 50),ROOT.TF1("pol1Fit", "pol1", 0, 50)]
	fitNames = ["pol2Fit", "expoFit", "pol1Fit"]
	#factorialFit = ROOT.TF1("factorialFit","[0]*Factorial(x)", 0, 50)
	#factorialFit.SetParameter(0, 1);

	for counter1 in range(6):
		graphs[counter1].GetXaxis().SetTitle("Number of Bunches")	
		graphs[counter1].GetYaxis().SetTitle("Rate")
		graphs[counter1].GetXaxis().SetLimits(0, 51)
		for counter2 in range(3):
			c = TCanvas("C", "Canvas")
			c.SetGrid()
			graphs[counter1].Fit(fits[counter2], "RQ")
			#graphs[counter1].Fit(factorialFit, "RQ")
			print graphs[counter1].GetName() + "_" + fitNames[counter2] + ", " + str(fits[counter2].Eval(2000)) + ", " + str(fits[counter2].Eval(1000)) + ", " + str(fits[counter2].GetChisquare())
			graphs[counter1].Draw("ALP")
			c.Modified()
			c.Update()			
			img = ROOT.TImage.Create()
			img.FromPad(c)
			img.WriteImage(graphs[counter1].GetName() + "_" + fitNames[counter2] + ".png")
			del c, img

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
	lessData = [[0,0,0,0,0,0] for x in range(5)]
	abvData = [[0,0,0,0,0,0] for x in range(5)]

	for bunchType in files:
		print "Analysing Bunch type -> " + str(bunchtypes[count])
		nums = [0, 0, 0, 0, 0, 0]
		lessValues = [0,0,0,0,0] # totPions, totTracks, true, partrue, fake
		abvValues = [0,0,0,0,0]
		fileCount = 1
		for fileName in bunchType:
			print "Processing File num -> " + str(fileCount)
			fileCount += 1
			lessValuesTemp, abvValuesTemp = track_validation(fileName)
			for counter in range(5):
				lessValues[counter] += lessValuesTemp[counter]
				abvValues[counter] += abvValuesTemp[counter]

		print("lessValues -> " + str(lessValues[0]) + ", " + str(lessValues[1]) + ", " + str(lessValues[2]) + ", " + str(lessValues[3]) + ", " + str(lessValues[4])) 
		print("abvValues -> " + str(abvValues[0]) + ", " + str(abvValues[1]) + ", " + str(abvValues[2]) + ", " + str(abvValues[3]) + ", " + str(abvValues[4])) 
		print("less(E,F) -> " + str(lessValues[2] / lessValues[0]) + ", " + str(lessValues[4] / lessValues[1]) + ", Abv(E,F) -> " + str(abvValues[2] / abvValues[0]) + ", " + str(abvValues[4] / abvValues[1]))
		lessData[count][0] = lessValues[2] / lessValues[0] #Efficiency
		lessData[count][1] = math.sqrt(((lessValues[2] / lessValues[0])*(1-(lessValues[2] / lessValues[0])))/(lessValues[0]))
		lessData[count][2] = lessValues[3] / lessValues[0] #parEfficiency
		lessData[count][3] = math.sqrt(((lessValues[3] / lessValues[0])*(1-(lessValues[3] / lessValues[0])))/(lessValues[0]))
		lessData[count][4] = lessValues[4] / lessValues[1] #Fake Rate
		lessData[count][5] = math.sqrt(((lessValues[4] / lessValues[1])*(1-(lessValues[4] / lessValues[1])))/(lessValues[1]))
		abvData[count][0] = abvValues[2] / abvValues[0]
		abvData[count][1] = math.sqrt(((abvValues[2] / abvValues[0])*(1-(abvValues[2] / abvValues[0])))/(abvValues[0]))
		abvData[count][2] = abvValues[3] / abvValues[0]
		abvData[count][3] = math.sqrt(((abvValues[3] / abvValues[0])*(1-(abvValues[3] / abvValues[0])))/(abvValues[0]))
		abvData[count][4] = abvValues[4] / abvValues[1]
		abvData[count][5] = math.sqrt(((abvValues[4] / abvValues[1])*(1-(abvValues[4] / abvValues[1])))/(abvValues[1]))

		count += 1

	draw_graphs(lessData, abvData, bunchtypes)

if __name__=='__main__':
	main()