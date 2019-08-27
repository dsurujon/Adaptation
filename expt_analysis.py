## DS 08/22/19
# Analyze filtered WGS results from multiple experiments together

import os
from optparse import OptionParser
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.style.use('bmh')
mycolors = mpl.cm.get_cmap("gist_rainbow", 25)
plt.switch_backend('Agg')

## Use example 
# python expt_analysis.py -i ExptSheet -o outputdir

options = OptionParser(usage='%prog -i ExptSheet -o outputdir', 
		description = "Specify Experiment Sheet to use, and the output directory to save new files in")
		

options.add_option("-i", "--exptsheet", dest="exptsheet", help="experiment sheet")
options.add_option("-o", "--outputdir", dest="outputdir", help="output directory")

annot = pd.read_csv("Annotation_3Strains_Cleaned.csv")
def make_mut_summary(expt_df):
	alldata = pd.DataFrame()
	allsummary = pd.DataFrame()
	n_expts = expt_df.shape[0]
	for i in range(n_expts):
		thistable = pd.read_csv(expt_df['File'][i])
		expt = expt_df['Name'][i]
		mergeby = expt_df['MergeBy'][i]
		# merge with annotation table
		annot_sub = annot[[mergeby, "Tag1", "Category1", "No..of.Categories", "Location.Tag", "Gene.Name", "Gene.Description"]]
		annot_sub.rename(columns={mergeby: "LocusTag"}, inplace=True)
		thistable_annot = thistable.merge(annot_sub, how="left", on="LocusTag")
		#add a column with the experiment name
		thistable_annot['Experiment']=expt
		# add a dummy column
		thistable_annot['AG']=1
		#concat to alldata
		alldata = pd.concat([alldata, thistable_annot], axis=0)
		
		total = thistable.shape[0]
		noncoding = sum(pd.isnull(thistable.LocusTag))
		coding = total-noncoding
		unique = thistable.LocusTag.nunique()
		
		snp = sum(thistable.Type=="SNP")
		sub = sum(thistable.Type=="SUB")
		ins = sum(thistable.Type=="INS")
		dlt = sum(thistable.Type=="DEL")
		
		summarydf = pd.DataFrame({'Experiment':[expt],
								  'Total':[total],
								 'Coding':[coding],
								 'Noncoding':[noncoding],
								 'Unique_LT':[unique],
								 'SNP': [snp],
								 'SUB': [sub],
								 'INS': [ins],
								 'DEL': [dlt]})
		allsummary = pd.concat([allsummary, summarydf], axis=0)
	return(allsummary, alldata)

def main():
	#read input args
	opts, args = options.parse_args()
	exptsheet = opts.exptsheet
	outdir = opts.outputdir
	
	# if output directory doesn't exist, make one
	if os.path.exists(outdir)==False:
		os.makedirs(outdir)
		
	# read the experiment sheet
	exptsheet_df = pd.read_csv(exptsheet)
	
	mut_summary, mut_annot_all = make_mut_summary(exptsheet_df)
	#reorder columns
	mut_summary = mut_summary[['Experiment', 'Noncoding', 'Coding', 'INS', 'DEL', 'SUB', 'SNP', 'Total', 'Unique_LT']]
	# write summary to file
	mut_summary_file = os.path.join(outdir, "Mutation_summary.csv")
	mut_summary.to_csv(mut_summary_file)
	# write annotated mutations file
	mut_annot_file = os.path.join(outdir, "Mutation_annotation_all.csv")
	mut_annot_all.to_csv(mut_annot_file)
	
	
	# plot the types of mutations (SNP/INS/DEL/SUB)
	plt.clf()
	ax = mut_summary[['SNP','INS','DEL','SUB']].iloc[mut_summary.Total.argsort()].plot(kind="bar", stacked=True, colormap=mycolors)
	ax.set_xticklabels(mut_summary[['Experiment']].iloc[mut_summary.Total.argsort(),0])
	ax.set_ylabel('Number of Mutations')
	plt.tight_layout()
	mut_type_plot_file = os.path.join(outdir, 'Mutation_type.svg')
	plt.savefig(mut_type_plot_file)
	
	# plot the locations of mutations (Coding/noncoding)
	plt.clf()
	ax = mut_summary[['Coding','Noncoding']].iloc[mut_summary.Total.argsort()].plot(kind="bar", stacked=True, colormap=mycolors)
	ax.set_xticklabels(mut_summary[['Experiment']].iloc[mut_summary.Total.argsort(),0])
	ax.set_ylabel('Number of Mutations')
	plt.tight_layout()
	mut_type_plot_file = os.path.join(outdir, 'Mutation_coding.svg')
	plt.savefig(mut_type_plot_file)
	
	## plot the Tags  - counts are MUTATIONS
	plt.clf()
	aglist_grouped = mut_annot_all.groupby(
		by=['Experiment', 'Tag1'],
		as_index=False
		).sum()
	aglist_grouped_wide = aglist_grouped.pivot(index='Experiment', columns='Tag1', values='AG')
	ax = aglist_grouped_wide.plot(kind="bar", stacked=True, colormap=mycolors)
	ax.set_ylabel('Number of Mutations')
	plt.tight_layout()
	AG_TAG_plot_file = os.path.join(outdir, 'AG_TAG_mutations.svg')
	plt.savefig(AG_TAG_plot_file)
	
	## plot the Categories - counts of MUTATIONS
	plt.clf()
	aglist_grouped = mut_annot_all.groupby(
		by=['Experiment', 'Category1'],
		as_index=False
		).sum()
	aglist_grouped_wide = aglist_grouped.pivot(index='Experiment', columns='Category1', values='AG')
	ax = aglist_grouped_wide.plot(kind="bar", stacked=True, colormap=mycolors)
	ax.set_ylabel('Number of Mutations')
	plt.tight_layout()
	AG_CAT_plot_file = os.path.join(outdir, 'AG_CATEGORY_mutations.svg')
	plt.savefig(AG_CAT_plot_file)
	
	
	aglist = mut_annot_all.drop_duplicates(['Experiment', 'LocusTag', 'Category1', 'Tag1'])
	# remove non-locus tags (intergenics) 
	aglist_genic = aglist.dropna(subset=['LocusTag'])
	## plot the Tags  - counts are GENES
	plt.clf()
	aglist_grouped = aglist_genic.groupby(
		by=['Experiment', 'Tag1'],
		as_index=False
		).sum()
	aglist_grouped_wide = aglist_grouped.pivot(index='Experiment', columns='Tag1', values='AG')
	ax = aglist_grouped_wide.plot(kind="bar", stacked=True, colormap=mycolors)
	ax.set_ylabel('Number of Genes')
	plt.tight_layout()
	AG_TAG_plot_file = os.path.join(outdir, 'AG_TAG_genes.svg')
	plt.savefig(AG_TAG_plot_file)
	
	## plot the Categories - counts of GENES
	plt.clf()
	aglist_grouped = aglist_genic.groupby(
		by=['Experiment', 'Category1'],
		as_index=False
		).sum()
	aglist_grouped_wide = aglist_grouped.pivot(index='Experiment', columns='Category1', values='AG')
	ax = aglist_grouped_wide.plot(kind="bar", stacked=True, colormap=mycolors)
	ax.set_ylabel('Number of Genes')
	plt.tight_layout()
	AG_CAT_plot_file = os.path.join(outdir, 'AG_CATEGORY_genes.svg')
	plt.savefig(AG_CAT_plot_file)

	

if __name__ == '__main__':
	main()