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

def getHitPositions(event):
	# Get the Position data for all hits in event, calculate r and theta also,
	# Add these to dict (posData).
	barrelHits = event.getCollection("TrackerHits")
	hitCount = 0
	posData = defaultdict(list)
	for hit in barrelHits:
		position = hit.getPosition() # Gets DIGITISED position of hit, in (mm).
		if math.sqrt(pow(position[0],2) + pow(position[1],2)) > 150:
			hitCount += 1
			posData[hitCount].append(position[0]) # x
			posData[hitCount].append(position[1]) # y
			posData[hitCount].append(position[2]) # z
			posData[hitCount].append(math.sqrt(pow(position[0],2) + pow(position[1],2))) # r
			posData[hitCount].append(math.atan(position[1] / position[0])) # theta
	return posData

def minFun(val1, val2, f, par, t):
	np = gr.GetN()
	f = 0
	i = 0
	x = gr.GetX()
	y = gr.GetY()
	while i < np:
		u = x[i] - par[0]
		v = y[i] - par[1]
		dr = par[2] - math.sqrt(pow(u,2)+pow(v,2))
		f += pow(dr,2)
		i += 1

def xyRes(posData, rmlayer):
	if len(posData) == 5:
		c = TCanvas("c", "BOOM!!!!")
		gr = ROOT.TGraph(4)
	for hit in posData:
		if layer != rmlayer:
			gr.SetPoint(layer, posData[layer][0], posData[layer][1])
	ROOT.TVirtualFitter.SetDefaultFitter("Minuit")
	fitter = ROOT.TVirtualFitter.Fitter(0, 3)
	fitter.SetFCN(minFun)
	fitter.SetParameter(0, "x0", 0, 0.1, 0, 0)
	fitter.SetParameter(1, "y0", 0, 0.1, 0, 0)
	fitter.SetParameter(2, "R", 1, 0.1, 0, 0)
	arglist = []
	fitter.ExecuteCommand("MIGRAD", arglist, 0)
	arc = ROOT.TArc(fitter.GetParameter(0),fitter.GetParameter(1),
		fitter.GetParameter(2))
	gr.Draw("p")
	arc.Draw()

	arc = ROOT.TArc(gMinuit.GetParameter(0),gMinuit.GetParameter(1),
		gMinuit.GetParameter(2))
	gr.Draw("p")
	arc.Draw()

def parparRes(posData, layer, par1, par2):
	#c = TCanvas("c", "boom")
	layers = [215.075, 465.075, 715.075, 965.075, 1215.075]
	rad = []
	z = []
	resHit = []
	residuals = []
	hitCount = 0
	if len(posData) < 8:
		for hit in posData:
			#hitCount += 1
			if posData[hit][3] < (layers[layer]-100) or posData[hit][3] > (layers[layer]+100):
				rad.append(posData[hit][par1])
				z.append(posData[hit][par2])
				hitCount += 1
			if posData[hit][3] > (layers[layer]-100) and posData[hit][3] < (layers[layer]+100):
				resHit.append(hit)

		radArray = array.array("f", rad)
		zArray = array.array("f", z)
		szGraph = ROOT.TGraph(hitCount, zArray, radArray)
		szGraph.Fit("pol1","Q")
		#szGraph.Draw()
		#raw_input()
		for res in resHit:
			zActual = posData[res][par2]
			zFit = (posData[res][par1] - szGraph.GetFunction("pol1").GetParameter(0)) / szGraph.GetFunction("pol1").GetParameter(1)
			residual = zFit - zActual
			if fabs(residual) < 1: 
				residuals.append(residual)

		szGraph.Delete()
		del(zArray)
		del(radArray)

	return residuals

def res_graph(residuals, name, bins, width):
	hist = ROOT.TH1F(name, name, bins, -width, width)
	for res in residuals:
		hist.Fill(res)
	hist.Fit("gaus")
	# print hist.GetFunction("gaus").GetParameter(2)
	# print hist.GetRMS()
	return float(hist.GetFunction("gaus").GetParameter(2)), float(hist.GetFunction("gaus").GetParError(2))

def myFun(x):
	return x / math.sqrt(12)

def resSizeGraph(name, res, resError, size):
	line = ROOT.TGraph()
	x = [0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20]
	lineValues = [0.0231, 0.0289, 0.0346, 0.0404, 0.0462, 0.0520, 0.0577]
	graph = ROOT.TGraphErrors()
	count = len(res)
	for index in range(count):
		residual = float(res[index])
		residualError = float(resError[index])
		pixelSize = float(size[index])
		graph.SetPoint(index, pixelSize, residual)
		graph.SetPointError(index, 0, residualError)
		line.SetPoint(index, x[index], lineValues[index])

	graph.GetXaxis().SetTitle("Square Pixel Size[mm]")
	graph.GetYaxis().SetTitle("#sigma(Residual Distribution)[mm]")
	graph.SetMinimum(0.02)
	graph.SetMaximum(0.1)
	graph.GetXaxis().SetLimits(0.075, 0.205)
	c = TCanvas("C", name)
	c.SetGrid()
	graph.Draw("AP")
	line.Draw("SAME LP")
	c.Modified()
	c.Update()
	raw_input()
	img = ROOT.TImage.Create()
	img.FromPad(c)
	img.WriteImage(name + ".jpeg")
	del c, img

def main():
	# main
	args = commonFunctions.simple_parse_args()

	if not args:
		print 'Invalid Arguments'
		sys.exit(1)

	fileList = commonFunctions.input_files(args.inputDir, ".slcio")

	fileCount = 0
	residuals = []
	resErrors = []
	sizes = []
	for fileName in fileList:
		residualsRZ = []
		residualsXY = []
		fileCount += 1
		print "File -> " + str(fileCount)

		outputName, eventCount, detector, energy = commonFunctions.get_fileData(fileName)

		reader = IOIMPL.LCFactory.getInstance().createLCReader()
		reader.open(fileName)

		for event in reader:
			posData = getHitPositions(event)
			residualsRZ.extend(parparRes(posData, 3, 3, 2))
			residualsXY.extend(parparRes(posData, 3, 1, 0))


		sigma, sigmaError = res_graph(residualsRZ, "res", 50, 0.1)
		size = fileName[73:-21]
		size = "0." + size
		print size
		residuals.append(sigma)
		resErrors.append(sigmaError)
		sizes.append(size)

		#commonFunctions.plot_simple_histogram(fileName + "_RZ", residualsRZ, 50, -0.05, 0.05, True, 0, 0.015)
		#commonFunctions.plot_simple_histogram(fileName + "_XY", residualsXY, 50, -1, 1, False, 10, 10)
	resSizeGraph("Residual vs Pixel Size", residuals, resErrors, sizes)
	print "[DONE]"

if __name__=='__main__':
	main()