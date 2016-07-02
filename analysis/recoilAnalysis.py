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

def get_event_type(event):
	MCCollection = event.getCollection("MCParticlesSkimmed")
	h_count = z_count = 0
	for particle in MCCollection:
		if particle.getPDG() == 23: z_count += 1
		if particle.getPDG() == 25: h_count += 1
	if h_count == 1 and z_count == 1:
		return 1 # denotes a ZH event!
	elif h_count == 0 and z_count == 2:
		return 2 # denotes a ZZ event!
	else: return 0
	# If no return i.e false shall not use the event. Not easily extended within this script!!!!

def COM_energies(event):
	MCCollection = event.getCollection("MCParticlesSkimmed")
	total_energy = reduced_energy = count = 0
	for particle in MCCollection:
		if particle.getParents().size() == 0:
			total_energy += particle.getEnergy()
			reduced_energy += particle.getDaughters()[0].getEnergy()
			count += 1
	if count == 2:
		return total_energy, reduced_energy

def PFO_processor(event):
	muPlus, muMinus = False, False
	plusList, minusList = [], []

	PFOCollection = event.getCollection("PandoraPFOCollection")
	for particle in PFOCollection:
		if particle.getType() == 13:
			muPlus = True
			plusList.append(particle)
		if particle.getType() == -13:
			muMinus = True
			minusList.append(particle)

	parentMom = parentMass = 0

	if not muMinus or not muPlus:
		return False, parentMass, parentMom

	massDiff = 10	
	rParams = []

	for plusMuon in plusList:
		for minusMuon in minusList:
			parentEnergy = plusMuon.getEnergy() + minusMuon.getEnergy()
			muPlusMom = plusMuon.getMomentum()
			muMinusMom = minusMuon.getMomentum()
			# muPlusMomT = math.sqrt(pow(muPlusMom[0],2)+pow(muPlusMom[1],2))
			# muMinusMomT = math.sqrt(pow(muMinusMom[0],2)+pow(muMinusMom[1],2))
			parentMom = math.sqrt(pow((muPlusMom[0]+muMinusMom[0]),2)+pow((muPlusMom[1]+muMinusMom[1]),2)+pow((muPlusMom[2]+muMinusMom[2]),2))
			parentMass = math.sqrt(pow(parentEnergy,2)-pow(parentMom,2))
			if math.fabs(parentMass - 91.1876) < massDiff:
				massDiff = math.fabs(parentMass - 91.1876)
				rParams = [parentMass, parentMom]


	if len(rParams) > 0:
		# print str(rParams[0]) + ", " + str(rParams[1]) + ", " +str(rParams[4])
		return True, rParams[0], rParams[1]

	else:
		return False, parentMass, parentMom

def event_processor(event):
	event_type = get_event_type(event)
	if event_type is not 0:
		total_energy, reduced_energy = COM_energies(event)
		check, parentMass, parentMom = PFO_processor(event)
		if check:
			total_recoil_mass = get_recoil_mass(total_energy, parentMass, parentMom)
			reduced_recoil_mass = get_recoil_mass(reduced_energy, parentMass, parentMom)
			return True, event_type, total_recoil_mass, reduced_recoil_mass, reduced_energy

	return False, event_type, 0, 0, 0

def get_recoil_mass(energy, parentMass, parentMom):
	parentEnergy = math.sqrt(pow(parentMass, 2) + pow(parentMom, 2))
	recoilEnergy = energy - parentEnergy
	recoilMass = math.sqrt(fabs(pow(recoilEnergy, 2) - pow(parentMom, 2)))
	# print "----Recoil Mass = " + str(recoilMass) 
	return recoilMass

def combined_histogram(name, zh_data, zz_data, ratio, bins, minX, maxX):
	print "Plotting-> " + name
	zh_hist = ROOT.TH1F("zh_hist", "zh_hist", bins, minX, maxX)
	zz_hist = ROOT.TH1F("zz_hist", "zz_hist", bins, minX, maxX)
	for zh in zh_data:
		zh_hist.Fill(zh)
	print len(zh_data)
	for zz in zz_data:
		zz_hist.Fill(zz)
	print len(zz_data)

	zz_hist.Scale(ratio)
	# zz_hist.Scale(60)
	zh_hist.Add(zz_hist)
	c = TCanvas("c", name)

	zh_hist.Fit("gaus","","",105,145)
	
	zh_hist.Draw()
	zh_hist.Print()

	c.Modified()
	c.Update()

	img = ROOT.TImage.Create()
	img.FromPad(c)
	img.WriteImage(name + ".png")
 
def main():
	# main
	args = commonFunctions.simple_parse_args()
	if not args:
		print 'Invalid Arguments'
		sys.exit(1)

	fileList = commonFunctions.input_files(args.inputDir, ".slcio")

	data = defaultdict(list)
	zh_count = zz_count = 0 
	for fileName in fileList:
		outputName, eventCount, detector, energy = commonFunctions.get_fileData(fileName)
		
		reader = IOIMPL.LCFactory.getInstance().createLCReader()
		reader.open(fileName)

		muonCount = 0
		for event in reader:
			check, event_type, total_recoil_mass, reduced_recoil_mass, reduced_energy = event_processor(event)
			if check:
				muonCount += 1
				if event_type == 1: zh_count += 1
				if event_type == 2: zz_count += 1

				data["totalMass " + str(event_type)].append(total_recoil_mass)
				data["reducedMass " + str(event_type)].append(reduced_recoil_mass)
				data["reducedEnergy " + str(event_type)].append(reduced_energy)

		print ", Total Events-> " + str(eventCount) + ", DiMuon Events-> " + str(muonCount) 

	if zz_count > 0:
		ratio = zh_count / zz_count
	else: ratio =0

	for key in data:
		# Makes plots of all the lists in the dictionary data[] with appropriae names.
		if key.split()[0] == "totalMass" or key.split()[0] == "reducedMass":
			if key.split()[1] == "1":
				commonFunctions.plot_simple_histogram(args.inputDir + key, data[key], args.numberOfBins, 50, 250, True, 120, 20)
			if key.split()[1] == "2":
				commonFunctions.plot_simple_histogram(args.inputDir + key, data[key], args.numberOfBins, 50, 250, True, 91, 30)

		# elif key.split()[0] == "reducedEnergy":
			# commonFunctions.plot_simple_histogram(args.inputDir + key, data[key], args.numberOfBins, 350, 500, False)
			'''
		for key2 in data:

			# Makes 2D plots of the different types of energies
			if key.split()[0] == "reducedMass" and key2.split()[0] == "reducedEnergy" and key.split()[1] == key2.split()[1]:
				commonFunctions.plot_2d_histogram(args.inputDir + "reduced_2d_" + key.split()[1], data[key], data[key2], args.numberOfBins, 50, 200, 350, 500, False)
			elif key.split()[0] == "totalMass" and key2.split()[0] == "reducedEnergy" and key.split()[1] == key2.split()[1]:
				commonFunctions.plot_2d_histogram(args.inputDir + "total_2d_" + key.split()[1], data[key], data[key2], args.numberOfBins, 50, 200, 350, 500, False)
			elif key.split()[0] == "totalMass" and key2.split()[0] == "totalMass" and key.split()[1] < key2.split()[1]:
				combined_histogram(args.inputDir + "totalMass_comb", data[key], data[key2], ratio, args.numberOfBins, 0, 250)
			elif key.split()[0] == "reducedMass" and key2.split()[0] == "reducedMass" and key.split()[1] < key2.split()[1]:
				combined_histogram(args.inputDir + "reducedMass_comb", data[key], data[key2], ratio, args.numberOfBins, 0, 250)
'''
	reader.close()
	
	print "[DONE]"


if __name__=='__main__':
	main()