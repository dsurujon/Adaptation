## DS 08/22/19
# Copy aerobio/breseq output gd files into a new directory
# and rename them from "output.gd" to "[sampleame].gd" 
# where [samplename] is the name of the parent directory 
# of the output.gd file under ExpOut

import os
import shutil
from optparse import OptionParser

## Use example 
# python copy_gd.py -e expID -o outputdir

options = OptionParser(usage='%prog -e expID -o outputdir',
		description = "Specify the aerobio experiment ID (where the gd files are currently) and the directory you want to copy them to. The output.gd files will automatically be renamed in the new directory as their sample name.")
		
options.add_option("-e", "--expID", dest="expID", help="aerobio experiment ID")
options.add_option("-o", "--outputdir", dest="outputdir", help="output directory")


def main():
	#read input args
	opts, args = options.parse_args()
	expID = opts.expID
	outdir = opts.outputdir
	
	# if output directory doesn't exist, make one
	if os.path.exists(outdir)==False:
		os.makedirs(outdir)
	
	# list all subdirectories under expID/Out
	aerobio_dir = "/store/data/ExpOut/"+expID+"/Out/"
	samples = os.listdir(aerobio_dir)
	
	# for each sample, copy the gd file to outdir, and rename it
	for samplename in samples:
		if samplename!="Bams":
			print("Transferring" , samplename)
			try:
				source_gd = os.path.join(aerobio_dir, samplename, "output", "output.gd")
				shutil.copy(source_gd, outdir)
				oldgdname = os.path.join(outdir, "output.gd")
				newgdname = os.path.join(outdir, samplename+".gd")
				os.rename(oldgdname,newgdname)
			except IOError:
				print("Error transferring file", samplename)

if __name__ == '__main__':
	main()