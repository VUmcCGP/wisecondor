#!/usr/bin/env python
# Daphne van Beek

import matplotlib
from pyPdf import PdfFileWriter, PdfFileReader
import sys
import os
import glob
import re
from pylab import *
from reportlab.lib import colors as ncolors
from reportlab.lib.pagesizes import letter, inch
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Image, Spacer, PageBreak
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.pdfgen import canvas

#########################
# Helper Definitions	#
#########################
def natural_key(string_):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]
    
def append_pdf(input,output):
	[output.addPage(input.getPage(page_num)) for page_num in range(input.numPages)]

#########################################
# Start processing each input directory	#
#########################################
usage = "Usage: python wc_makeReport.py [inputdir 1] [inputdir 2]...[inputdir n] [outputdir]\nInput directories should contain per sample at least:\n-Wisecondor [samplename].txt and [samplename].pdf output files\n-[samplename].sort.bam.flagstat file resulting from wc_qc.sh"
nrOfDirs = len(sys.argv) - 2
outputDir = sys.argv[-1]
if len(sys.argv) < 3:
	print usage
	sys.exit()
print outputDir
if not os.path.exists(outputDir):
	os.makedirs(outputDir)

#########################
# Set (style) defaults	#
#########################
styles = getSampleStyleSheet()
styleHeading = styles['Heading2']
styleHeading.alignment = 1
styleNormal = styles['Heading4']
styleNormal.alignment = 1
styleNormal2 = styles["Normal"]
styleNormal2.alignment = 2
elements = []

Ncols = 1
Nrows = 1
edge = 2.5
figheight = edge * Nrows   # inches
figwidth = edge * Ncols    # inches
rcParams['font.size'] = 6.0
rcParams['axes.titlesize'] = 12.0
rcParams['xtick.labelsize'] = 6.0
rcParams['ytick.labelsize'] = 6.0
rcParams['legend.fontsize'] = 5.0
plotheight = figwidth/Ncols
H = 1.0	/ Nrows	#plotheight/figheight
W = 2.0 / Ncols
margin = 0.25
widthplot = W*(1-2*margin)
heightplot = H*(1-2*margin)

piePos = 1
figno = 0

for i in range(1,(nrOfDirs + 1)):
	path = sys.argv[i]
	print os.path.basename(os.path.normpath(path))
	print 'Samplename\tTotal\tMapped\tFiltered\tMQ0\tUnmapped\tLeft after filtering'

	
	#################################################
	# Create report for each sample in directory	#
	#################################################
	for files in glob.glob(path + '*.txt'):
		resultLines = []
		name = os.path.basename(os.path.normpath(files))
		name = name.split(".")[0]
		pathSample = files.split(".")[0:-1]
		pathSample = ".".join(pathSample)
		outputPDF = SimpleDocTemplate(outputDir + "/" + name + '_report.pdf', pagesize=letter)

		#########################################
		# Read information from input files		#
		#########################################
		flagstat = open(pathSample + ".sort.bam.flagstat", 'r')
		flTotal = float(flagstat.next().split(" ")[0])
		flDups = int(flagstat.next().split(" ")[0])
		flMapped = float(flagstat.next().split(" ")[0])
		for i in range(1,10):
			flagstat.next()
		flMQ0 = int(flagstat.next().split(" ")[-1])
		flagstat.next()
		flRETRO = float(flagstat.next().rstrip())
		flagstat.close()

		'''
		configFile = open(glob.glob(pathSample + ".conf*")[0], 'r')
		rodaVersion = configFile.next().split()[2:]
		dateRun = configFile.next().split()[-1]
		configFile.close()
		'''	
			
		IN = open(files, 'r')
		IN.next(); IN.next()
		binsize = IN.next().rstrip().split("\t")[-1]
		for j in range(0,25):
			IN.next()
		reference = IN.next().rstrip().split("\t")[-1]
		for j in range(0,8):
			IN.next()
		windowsize = IN.next().rstrip().split("\t")[-1]
		for j in range(0,10):
			IN.next()
		averageAllowedDev = IN.next().rstrip().split(":")[1]
		averageAllowedDev = averageAllowedDev[1:6]
		
		binTest = []
		chrTest = False
		for line in IN:
			if str(line[0:10]) == "# Results:":
				IN.next()
				line = IN.next().rstrip()
				while line[0:21] != "# Script information:":
					if line[0:11] == "Single bin,":
						line = line[0:-1] + " (blue line):"
					elif line[0:9] == "Windowed,":
						line = "Smoothed/" + line[0:-1] + " (red line):"
					elif chrTest == True:
						linetest = line.split("\t")
						if len(linetest) > 2:
							if float(linetest[2]) < 4 and float(linetest[2]) > -4:
								line = ""
							else:
								line = "\tCall:\t" + linetest[1] + "\tZ-score:\t" + linetest[2]
					elif "+" in line:
						i = line.find("+")
						line = line[0:i] + "Dup." + line[i+1:]
					elif "-" in line:
						i = line.find("-")
						if line[i+1] == "\t":
							line = line[0:i] + "Del." + line[i+1:]
					elif "marked" in line:
						line = line[0:-6] + "of chromosome is indicated as duplicated/deleted"
					elif line[0:33] == "Chromosome wide, aneuploidy test:":
						line = line[0:-1] + " (based on red line):"
						chrTest = True						
					resultLines.append(line)
					line = IN.next().rstrip()
				break
		
		resultLines.append("Average allowed devation (AvgASD):")
		resultLines.append("\t" + str(averageAllowedDev) + " %")
		resultLines.append("Bin size:")
		resultLines.append("\t" + str(binsize) + " nt")
		resultLines.append("Window size:")
		resultLines.append("\t" + str(windowsize) + " bins")
		resultLines.append("Reference used:")
		resultLines.append("\t" + str(reference))
		resultLines.append("")
		
		#########################################
		# Create Mapped vs. Unmapped pie chart	#
		#########################################
		figno += 1
		figure(figno, figsize=(figwidth, figheight))
		flUnmapped		= 	flTotal - flMapped
		flFiltered = flMapped - flMQ0 - flRETRO
		percUnmapped	=	100 * flUnmapped / flTotal
		percMQ0			=	100 * flMQ0 / flTotal
		percMappedUsed	=	100 * (flMapped - flMQ0 - flFiltered) / flTotal
		percFiltered	=	100 * flFiltered / flTotal

		labels=['Mapped','RETRO Filter', 'MQ0','Unmapped']
		fracs = [percMappedUsed,percFiltered,percMQ0,percUnmapped]
		legend1 = [labels[i] for i in range(len(labels))]
		explode=(0.0, 0.0, 0.0, 0.0)
		colors=('#34C934','#FF9900','#FFFF00','#FF0000')
		#axes([left[piePos], bottom[piePos], widthplot, heightplot])
		patches, texts, autotexts = pie(fracs, labels=labels,colors=colors, autopct='%1.f%%', shadow=False) #,  explode=explode
		for ent in range(0,len(autotexts)):
			autotexts[ent].set_fontsize(6)
		title('Mapping statistics')
		legend(legend1,loc=(-0.1,-0.1))
		piePos += 1

		savefig(pathSample + '_pie.png', dpi=300)

		labels = []
		fracs = []
		legend1 = []
		#########################################
		# Print metrics statistics to tables	#
		#########################################
		#elements.append(Paragraph('RoDa Version ' + ' '.join(rodaVersion) + ' Quality report', styleHeading))
		elements.append(Paragraph('Quality report', styleHeading))
		elements.append(Paragraph('Sample ' + name, styleHeading))
		elements.append(Spacer(inch, 0.1*inch)	)
		now = datetime.datetime.now()
		#elements.append(Paragraph('Date analysis performed:\t' + str(dateRun), styleNormal2))
		elements.append(Paragraph('Date printing report:\t' + now.strftime("%d/%m/%y"), styleNormal2))
		elements.append(Spacer(inch, 0.25*inch)	)	
		elements.append(Paragraph('Mapping quality report', styleNormal))				
		
		[flTotal, flMapped, flMQ0, flUnmapped, flFiltered, flRETRO] = ('%d'%i if i == int(i) else '%s'%i for i in [flTotal, flMapped, flMQ0, flUnmapped, flFiltered, flRETRO])
		print name + '\t' + str(flTotal) + '\t' + str(flMapped) + '\t' + str(flFiltered) + '\t' + str(flMQ0) + '\t' + str(flUnmapped) + '\t' + str(flRETRO)

		data = [['Total reads', 					str(flTotal), 		str(100) + ' %'],
				['Mapped reads', 					str(flMapped),		str(round((float(flMapped) / float(flTotal) * 100),2)) + ' %'],
				['MQ0 reads', 						str(flMQ0),			str(round((float(flMQ0) / float(flTotal) * 100),2)) + ' %'],
				['Unmapped reads', 					str(flUnmapped), 	str(round((float(flUnmapped) / float(flTotal) * 100),2)) + ' %'],
				['Reads filtered by RETRO',			str(flFiltered),	str(round((float(flFiltered) / float(flTotal) * 100),2)) + ' %'],
				['Reads used by WISECONDOR',		str(flRETRO), 		str(round((float(flRETRO) / float(flTotal) * 100),2)) + ' %']]
			
		t = Table(data, 3 * [2.05 * inch], len(data) * [0.2 * inch])
		t.setStyle(TableStyle([	('INNERGRID', (0,0), (-1,-1), 0.25, ncolors.black),
								('BOX', (0,0), (-1,-1), 0.25, ncolors.black),
								('VALIGN', (0,0), (-1,-1), 'MIDDLE')]))
								
								
		elements.append(t)
		elements.append(Spacer(inch, 0.65*inch)	)	
		elements.append(Image(pathSample + '_pie.png', figwidth*inch, figheight*inch))

		rcParams['font.size'] = 9.0
		rcParams['xtick.labelsize'] = 6.0
		
		#################################
		# Print the Wisecondor Results	#
		#################################
		
		linePos = 700
		lineSpacing = 14
		origRowPos = 100
		c = canvas.Canvas(outputDir + "/" + name + "_results.pdf", pagesize=letter)
		cwidth, cheight = letter
		c.drawString(origRowPos,linePos,"RESULTATEN WISECONDOR")
		linePos = linePos - lineSpacing
		c.line(origRowPos, linePos, cwidth - origRowPos, linePos)
		for line in resultLines:
			rowPos = origRowPos
			tabSpacing = 50
			if line:
				if line == "Smoothed/Windowed, bin test (red line):" or line == "Average allowed devation (AvgASD):" or line == "Chromosome wide, aneuploidy test (based on red line):":
					c.line(rowPos, linePos, cwidth - rowPos, linePos)
					linePos = linePos - lineSpacing
				if line[-1] == ":":
					#linePos = linePos - lineSpacing
					c.drawString(rowPos,linePos,"")
					c.setFont("Helvetica-Bold",10)
					linePos = linePos - lineSpacing
					c.drawString(rowPos,linePos,line)
					linePos = linePos - lineSpacing
				else:
					linet = line.split("\t")
					c.setFont("Helvetica",10)
					#linePos = linePos - lineSpacing
					if len(linet) > 1:
						for tabs in linet:
							c.drawString(rowPos,linePos,tabs)
							rowPos = rowPos + tabSpacing
					else:					
						c.drawString(rowPos,linePos,line)
					linePos = linePos - lineSpacing
		
		c.line(origRowPos, 40, cwidth - origRowPos, 40)
		c.drawString(origRowPos, 30, "Bron: R. Straver et al. WISECONDOR: Detection of fetal aberrations from shallow sequencing") 
		c.drawString(origRowPos, 20, "maternal plasma based on a within-sample comparison scheme. 2013")

		c.save()
		
		
		#############################
		# Export the pdf document	#
		#############################
		outputPDF.build(elements)
		os.remove(pathSample + '_pie.png')
		
		#############################################
		# Concatenate Wisecondor plot and report	#
		#############################################
		#elements.append(Image(pathSample + ".zscore.pdf"))
		output = PdfFileWriter()
		append_pdf(PdfFileReader(file(outputDir + "/" + name + "_report.pdf","rb")),output)
		append_pdf(PdfFileReader(file(outputDir + "/" + name + "_results.pdf","rb")),output)
		append_pdf(PdfFileReader(file(pathSample + ".zscore.pdf","rb")),output)
 
		output.write(file(outputDir + "/" + name + ".pdf","wb"))
		os.remove(outputDir + "/" + name + "_results.pdf")
		os.remove(outputDir + "/" + name + "_report.pdf")
			
