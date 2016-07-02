# THIS SCRIPT USES THE NEW API
# author: jt12194@my.bristol.ac.uk based on script by jan.strube@cern.ch
from DIRAC.Core.Base import Script
Script.initialize()
from ILCDIRAC.Interfaces.API.DiracILC import DiracILC
from ILCDIRAC.Interfaces.API.NewInterface.UserJob import UserJob
from ILCDIRAC.Interfaces.API.NewInterface.Applications import SLIC, LCSIM, SLICPandora, Marlin, OverlayInput
from DIRAC.Resources.Catalog.FileCatalogClient import FileCatalogClient
import sys
import os.path
import os.path as path
import re
import argparse

def parse_args():
	parser = argparse.ArgumentParser(description='Handle the full SiD reconstruction chain')

	# The .stdhep input LFN, must be a file on the GRID.
	parser.add_argument('Input', help= 'LFN path to the input file')

	# The total number of arguments you want to run from this file. You can find the number of events in a .stdhep file by running,
	# "/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-08/lcio/v02-06/bin/lcio_event_counter file.stdhep" on an lxplus node.
	parser.add_argument('-r' ,'--events', help='Total number of events to use from input file', default = -1)

	parser.add_argument('-s','--split', help='The number of events to output into seperate output files (must be a value less than --runs). Default not split'
					  , default = -1)

	parser.add_argument('-c','--chain', nargs='+', type=int, help='What to run? In Order!!!\n1=Simulation\n2=AddOverlay\n3=digi/tracking\n4=slicPandora\n5=Vertexing/DST\n6=Flavortagging', required=True)

	parser.add_argument('-d', '--detector', help='Default is "sidloi3", give URL to .zip of new geometry (default=%(default)s)', default='sidloi3')

	parser.add_argument('--digiSteering', help='Path to digi steering file, If you want to include pixels etc...', default="steeringFiles/sid_dbd_prePandora_noOverlay.xml")

	parser.add_argument('-m', '--macFile', help='macro file to pass to SLIC (default=%(default)s)', default='steeringFiles/defaultILCCrossingAngle.mac')

	parser.add_argument('-b', '--bunches', help='Number of bunches to Overlay.', default = 1)

	parser.add_argument('--dontPromptMe', help='do not wait for user to confirm first job', action='store_true')

	return parser.parse_args()

def check_events_arguments(events, split_size):
	# Check the --runs and --split arguments to make sure they are compatible, if not exit...
	checkEvents = int(events) > 0 
	checkSplit = int(split_size) > 0

	if not checkEvents:
		print 'Error: The number of events must be specified with (--events) when calling script'
		sys.exit(1)

	elif checkEvents and checkSplit and int(split_size) > int(events):
		print 'Error: The total number of events (--runs) must be greater than the output size (--split)'
		sys.exit(1)

	elif checkEvents and not checkSplit:
		print 'Total number of events to run = ', events, '. NO SPLIT'

	elif checkEvents and checkSplit and int(events) > int(split_size):
		print 'Total number of events to run = ', events, '. SPLIT, into files of size = ', split_size, ' events.'

	else: 
		print "No Idea whats going on!!!!"
		sys.exit(1)

def input_output(input_lfn, detector, chain, digiSteering):
	lfn, extension = os.path.splitext(input_lfn)
	print 'Input LFN = ', input_lfn
	if 2 in chain and 3 not in chain:
		print "Error: Can only apply Overlay if following lcsim step (3) included"
		sys.exit(1)
	if 1 in chain and extension in ['.stdhep']:
		print 'Input LFN = ', input_lfn
	elif 1 in chain and extension not in ['.stdhep']:
		print 'Error: Need .stdhep file as input.'
		sys.exit(1)
	elif 2 in chain or 3 in chain or 4 in chain or 5 in chain or 6 in chain and extension in ['.slcio']:
		print 'Input LFN = ', input_lfn
	elif 2 in chain or 3 in chain or 4 in chain or 5 in chain or 6 in chain and extension not in ['.slcio']:
		print 'Error: Need .slcio file as input.'
		sys.exit(1)
	
	if len(chain) == 0: 
		print "Error: No processes to run."
		sys.exit(1)

	if digiSteering == "steeringFiles/sid_dbd_prePandora_noOverlay.xml":
		outputBase = path.basename(lfn) + '_' + detector + ".slcio"
		outputPath = detector
		repoName = "repo/" + path.basename(lfn) + '_' + detector + '_repository.cfg'
	if digiSteering != "steeringFiles/sid_dbd_prePandora_noOverlay.xml":
		if os.path.isfile(digiSteering):
			steeringBase, steeringExt = os.path.splitext(digiSteering)
			outputBase = path.basename(lfn) + '_' + detector + "_" + path.basename(steeringBase) + ".slcio"
			outputPath = detector + "_" + path.basename(steeringBase)
			repoName = "repo/" + path.basename(lfn) + '_' + detector + '_' + path.basename(steeringBase) + '_repository.cfg'
		else:
			print "Error: Steering File does not exist."
			sys.exit(1)

	outputPath = "overlayNew" # JUST PUT IN FOR THE OVERLAY STUDY!!! TAKE OUT!!!!!
	print "Output BaseName = ", outputBase
	print "Ouput Path = ", outputPath
	print "Repo Name = ", repoName
	return outputPath, outputBase, repoName

def setup_sandboxes(macFile): 
	# Need to have the pandoraSettings.xml on the grid or slicPandora fails for some reason. Also add Olivers lcsim driver.
	inputSandbox = ["LFN:/ilc/user/j/jtingey/pandoraSettings/pandoraSettings.xml"]
	inputSandbox.append("lcsimDrivers.jar")
	#inputSandbox.append("sidloi3FullPix.lcdd") # THIS NEED TO BE REMOVED AND FIXED!!!!
	outputSandbox = ['*.log', '*.xml'] # WOULD BE COOL TO HAVE THIS DEPENDANT ON INPUT	
	if macFile != 'defaultIlcCrossingAngle.mac':
		inputSandbox.append(macFile)
	return inputSandbox, outputSandbox

def alias_properties(detector):
	# Changes Alias Properties file, depending on the detector being used.
	aliasFile = open("steeringFiles/alias.properties", "w")
	aliasFile.truncate()
	line = detector + ': http://www.cern.ch/test-jtingey/' + detector + '.zip'
	print "Alias Properites line: " + line
	aliasFile.write(line)
	aliasFile.close()

def main():
	# Take the input arguments from the argument parser, and check they exist...
	args = parse_args()
	if not args:
		print 'Invalid Arguments'
		sys.exit(1)

	print args.chain[0]

	# softVersions = ["v3r0p3", "3.0-SNAPSHOT", "ILC_DBD", "0116"]
	softVersions = ["v3r0p3", "HEAD", "ILC_DBD", "0116"] # Working (recommended)
	# softVersions = ["v3r0p3", "2.5", "ILC_DBD", "0116"] # Working 
	# softVersions = ["v3r0p3", "HEAD", "ILC_DBD", "ILCSoft-01-17-07"] # Working
	# softVersions = ["v3r0p3", "HEAD", "ILCSoft-01-17-08", "0116"]

	check_events_arguments(args.events, args.split)
	detector = args.detector
	alias_properties(detector)
	outputPath, outputBase, repoName = input_output(args.Input, detector, args.chain, args.digiSteering)
	inputSandbox, outputSandbox = setup_sandboxes(args.macFile)

	dirac = DiracILC(True, repoName)

	# Prepares values for the job loop...
	if args.split < 0:
		nInputEvents = int(args.events)
		nOutputEvents = int(args.events)
	if args.split > 0:
		nInputEvents = int(args.events)
		nOutputEvents = int(args.split)

	# Loop that runs through the required number of jobs to be executed...
	for startEvent in range(0, nInputEvents, nOutputEvents):

################## Job Initialise ########################################		
		job = UserJob()
		job.setName(outputBase)
		job.setJobGroup('JobGroup')
		job.setInputSandbox(inputSandbox)
		fileNumber = startEvent/nOutputEvents
		print "Job ---> ", fileNumber

################## SLIC ##################################################
		if 1 in args.chain:
			slic = SLIC()
			slic.setVersion(softVersions[0])
			slic.setSteeringFile(args.macFile)
			# slic.setInputFile(lfn)
			slic.setOutputFile(outputBase.replace('.slcio', '_' + str(fileNumber) + '_sim.slcio'))
			slic.setDetectorModel(detector)
			slic.setNumberOfEvents(nOutputEvents)
			slic.setStartFrom(startEvent)
			#print slic.listAttributes()
			result = job.append(slic)
			if not result['OK']:
				print result['Message']
				sys.exit(2)

################## Overlay ###############################################
		if 2 in args.chain:
			'''
			#Add the gghad background overlay.
			gghad = OverlayInput()
			#gghad.setProdID(1767)
			gghad.setEnergy(500.0)
			gghad.setBXOverlay('args.bunches')
			gghad.setGGToHadInt( 4.1 )
			gghad.setNbSigEvtsPerJob(nOutputEvents)
			gghad.setMachine('ilc_dbd')
			gghad.setDetectorModel('sidloi3')
			gghad.setBkgEvtType('aa_lowpt')
			result = job.append( gghad )
			if not result['OK']:
				print result['Message']
				sys.exit(2)
			
			#Add the pair background overlay.
			pairs = OverlayInput()
			pairs.setProdID(2)
			pairs.setEnergy(500.0)
			pairs.setBXOverlay('args.bunches')
			pairs.setGGToHadInt(1.)
			pairs.setNbSigEvtsPerJob(nOutputEvents)
			pairs.setMachine('ilc_dbd')
			pairs.setDetectorModel('sidloi3')
			pairs.setBkgEvtType('eepairs')
			result = job.append( pairs )
			if not result['OK']:
				print result['Message']
				sys.exit(2)
			'''
			gghad = OverlayInput()
			gghad.setPathToFiles('/ilc/user/j/jstrube/gghadron_lowpt/sidloi3/')
			gghad.setBXOverlay(int(args.bunches))
			gghad.setGGToHadInt( 4.1 )
			gghad.setNbSigEvtsPerJob(nOutputEvents)
			gghad.setBkgEvtType('aa_lowpt')		
			result = job.append( gghad )
			if not result['OK']:
				print result['Message']
				sys.exit(2)
			'''
			pairs = OverlayInput()
			pairs.setPathToFiles('/ilc/user/j/jstrube/GuineaPig/sidloi3/')
			pairs.setBXOverlay(int(args.bunches))
			pairs.setGGToHadInt(1.)
			pairs.setBkgEvtType('eepairs')
			pairs.setNbSigEvtsPerJob(nOutputEvents)
			result = job.append( pairs )
			if not result['OK']:
				print result['Message']
				sys.exit(2)
			'''
			
################## lcsim (digitization and tracking) #####################
		if 3 in args.chain:
			lcsim = LCSIM()
			lcsim.setVersion(softVersions[1])
			lcsim.setSteeringFile(args.digiSteering) # Another version is included in /steeringFiles
			if 1 in args.chain:
				lcsim.getInputFromApp(slic)
			lcsim.setTrackingStrategy('steeringFiles/sidloi3_trackingStrategies_default.xml')
			lcsim.setAliasProperties('steeringFiles/alias.properties')
			lcsim.setDetectorModel(detector+".zip")
			#lcsim.setOutputFile(outputBase.replace('.slcio', '_' + str(fileNumber) + '_digiTracking.slcio'))
			lcsim.setOutputDstFile(outputBase.replace('.slcio', '_' + str(fileNumber) + '_DST.slcio')) #NEED TO CHANGE!!!
			lcsim.setNumberOfEvents(nOutputEvents)
			#print lcsim.listAttributes()
			result = job.append(lcsim)
			if not result['OK']:
				print result['Message']
				sys.exit(2)

################## slicPandora ###########################################
		if 4 in args.chain:
			slicPandora = SLICPandora()
			slicPandora.setVersion(softVersions[2])
			slicPandora.setDetectorModel(detector)
			slicPandora.getInputFromApp(lcsim)
			slicPandora.setOutputFile(outputBase.replace('.slcio', '_' + str(fileNumber) + '_pandora.slcio'))
			slicPandora.setPandoraSettings('pandoraSettings.xml')
			slicPandora.setNumberOfEvents(nOutputEvents)
			#print slicPandora.listAttributes()
			result = job.append(slicPandora)
			if not result['OK']:
				print result['Message']
				sys.exit(2)

################## Marlin, LCFIPlus Vertexing ############################
		if 5 in args.chain:
			vertexing = Marlin()
			vertexing.setVersion(softVersions[3])
			vertexing.setSteeringFile('steeringFiles/sid_dbd_vertexing.xml')
			vertexing.setGearFile('steeringFiles/' + detector + '.gear')
			vertexing.getInputFromApp(slicPandora)
			vertexing.setOutputFile(outputBase.replace('.slcio', '_' + str(fileNumber) + '_vertexing.slcio'))
			vertexing.setNumberOfEvents(nOutputEvents)
			#print vertexing.listAttributes()
			result = job.append(vertexing)
			if not result['OK']:
				print result['Message']
				sys.exit(2)
################## lcsim (DST production) ################################
			lcsimDst = LCSIM()
			lcsimDst.setVersion(softVersions[1])
			lcsimDst.setSteeringFile('steeringFiles/sid_dbd_postPandora.xml')
			lcsimDst.getInputFromApp(vertexing)
			lcsimDst.setNumberOfEvents(nOutputEvents)
			lcsimDst.setAliasProperties('steeringFiles/alias.properties')
			lcsimDst.setDetectorModel(detector+".zip")
			lcsimDst.setOutputRecFile(outputBase.replace('.slcio', '_' + str(fileNumber) + '_Rec.slcio'))
			lcsimDst.setOutputDstFile(outputBase.replace('.slcio', '_' + str(fileNumber) + '_DST.slcio'))
			#print lcsimDst.listAttributes()
			result = job.append(lcsimDst)
			if not result['OK']:
				print result['Message']
				sys.exit(2)

################## Marlin, LCFIPlus flavortag ############################
		if 6 in args.chain:
			flavortag = Marlin()
			flavortag.setVersion(softVersions[3])
			flavortag.setSteeringFile('steeringFiles/sid_dbd_flavortag.xml')
			flavortag.setGearFile('steeringFiles/' + detector + '.gear')
			flavortag.setInputFile(lcsimDstOutput)
			flavortag.setOutputFile(outputBase.replace('.slcio', '_' + '_flavortag.slcio'))
			flavortag.setNumberOfEvents(nOutputEvents)
			#print flavortag.listAttributes()
			result = job.append(flavortag)
			if not result['OK']:
				print result['Message']
				sys.exit(2)

################## Job Finalise ##########################################

		# List of banned sites that the job shall not be sent too. These are sites that jobs tend to fail on,
		# This list is likely to change.
		job.setBannedSites(['LCG.IN2P3-CC.fr', 'LCG.RAL-LCG2.uk', 'LCG.DESY-HH.de', 'LCG.DESYZN.de', 'LCG.KEK.jp',
							'OSG.PNNL.us','OSG.UConn.us','OSG.GridUNESP_CENTRAL.br','LCG.SCOTGRIDDURHAM.uk',
							'LCG.TECHNIONself.il','LCG.UKI-SOUTHGRID-RALPP.uk','OSG.FNAL_FERMIGRID.us','LCG.UKI-LT2-IC-HEP.uk'])

		job.setCPUTime(50000)
		job.setPlatform('x86_64-slc5-gcc43-opt')

		# Sets the output data file according to if -f is selcted, ships ouput to your /ilc/user/a/aPerson/
		# directory on the grid.
		outputLevel = max(args.chain)
		if outputLevel == 1:
			job.setOutputData(outputBase.replace('.slcio', '_' + str(fileNumber) + '_sim.slcio'), outputPath, 'CERN-SRM')
		if outputLevel == 3:
			#job.setOutputData(outputBase.replace('.slcio', '_' + str(fileNumber) + '_digiTracking.slcio'), outputPath, 'CERN-SRM')
			job.setOutputData(outputBase.replace('.slcio', '_' + str(fileNumber) + '_DST.slcio'), outputPath, 'CERN-SRM')
		if outputLevel == 4:
			job.setOutputData(outputBase.replace('.slcio', '_' + str(fileNumber) + '_pandora.slcio'), outputPath, 'CERN-SRM')
		if outputLevel == 5:
			job.setOutputData(outputBase.replace('.slcio', '_' + str(fileNumber) + '_DST.slcio'), outputPath, 'CERN-SRM')
		if outputLevel == 6:
			job.setOutputData(outputBase.replace('.slcio', '_' + str(fileNumber) + '_flavortag.slcio'), outputPath, 'CERN-SRM')

		job.setOutputSandbox(outputSandbox)
		job.setInputData(args.Input)

		if args.dontPromptMe:
			job.dontPromptMe()
		# Submits Job!!!
		job.submit()

	return 0;

if __name__=='__main__':
	main()