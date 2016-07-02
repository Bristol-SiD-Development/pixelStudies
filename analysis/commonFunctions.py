import os, sys, argparse, os.path
import ROOT
import array
import math

from pyLCIO import IOIMPL
from pyLCIO import UTIL
from pyLCIO import EVENT
from pyLCIO import IOIMPL

from ROOT import *

def input_files(inputDirectory, ext):
	# Checks the input files and returns list of those to process.
	inputFiles = []
	if os.path.isdir(inputDirectory):
		print '\nLooking in "' + inputDirectory + '" for input files...'
		fileNum = 0
		for fileName in os.listdir(inputDirectory):
			name, extension = os.path.splitext(fileName)
			if extension not in [ext]:
				print '----The file "' + fileName + '" is not valid, skipping!'

			else:
				print '----The File "' + fileName + '" is valid - adding!'
				inputFiles.append(os.path.join(inputDirectory, fileName))
				fileNum = fileNum + 1 

		if fileNum > 0:
			print "\nA total of " + str(fileNum) + " .slcio files shall be processed"
			return inputFiles

		else:
			print 'ERROR: No input files could be found!!!'
			sys.exit(1)
		
	else:
		print 'ERROR: Can only define either an input directory or an input file. Ensure they exist and are of the .slcio type!!!' 
		sys.exit(1)

def simple_parse_args():
	# Takes in arguments from the command line when script is called and returns them to main().
	currentDir = os.getcwd()
	parser = argparse.ArgumentParser(description='Processes .slcio file...'
						  ,epilog='In case of questions or problems, contact jt12194@my.bristol.ac.uk')

	parser.add_argument('inputDir', help='Input .slcio fileDirectory')

	parser.add_argument('-o', '--outputDirectory',
						help='Output directory for plots to be saved to, default = current directory.',
						default=currentDir)

	parser.add_argument('-bins', '--numberOfBins',
						help='The number of bins for the histograms.',
						default=100)

	return parser.parse_args()

def get_fileData(fileName):
	# Takes a file and prints/returns usefull info about it.
	inputName, inputExtension = os.path.splitext(fileName)
	outputName = os.path.basename(inputName)

	reader = IOIMPL.LCFactory.getInstance().createLCReader()
	reader.open(fileName)
	eventCounter = 0
	for event in reader:
		if eventCounter == 0:
			detector = event.getDetectorName()
			MCParticle = event.getCollection("MCParticlesSkimmed")
			energy = 0
			for particle in MCParticle:
				if particle.getParents().size() == 0:
					energy += particle.getEnergy()	
		eventCounter += 1

	reader.close()

	print "----Processing -> " + outputName
	print "NO.Events-> " + str(eventCounter) + ", Detector-> " + detector + ", EventEnergy-> " + str(energy)
	return outputName, eventCounter, detector, energy

def plot_simple_histogram(name, data, bins, minX, maxX, fit, centreFit, rangeFit):
	print "Plotting-> " + name 
	# Plots histogram and fits anything according to arguments.
	hist = ROOT.TH1F(name, name, bins, minX, maxX)
	c = TCanvas("c", name)
	for dataPoint in data:
		hist.Fill(dataPoint)
	
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


def plot_2d_histogram(name, data1, data2, bins, minX, maxX, minY, maxY, fit):
	print "Plotting-> " + name
	hist = ROOT.TH2F(name, name, bins, minX, maxX, bins, minY, maxY)
	c = TCanvas("c", name)
	for i in range(0, len(data1)):
		hist.Fill(data1[i], data2[i])

	hist.Draw("cont")
	hist.Print()

	if fit: print "Currently fitting not implemented"

	img = ROOT.TImage.Create()
	img.FromPad(c)
	img.WriteImage(name + ".png")

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


