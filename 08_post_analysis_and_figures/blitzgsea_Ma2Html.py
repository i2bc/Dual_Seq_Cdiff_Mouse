import blitzgsea as blitz
import urllib.request
import pandas as pd
from matplotlib import pyplot as plt

# read signature as pandas dataframe
sig = pd.read_csv("stat_CDS_bg.tsv")
signature=sig.dropna()  # suppress NA value
# use enrichr submodule to retrieve gene set library
library = blitz.enrichr.read_gmt("Ma2Html_for_blitzgsea_analysis", verbose=True)
# run enrichment analysis
result = blitz.gsea(signature, library, verbose = True, plotting = True, signature_cache = False, seed = -1, permutations = 10000)

# saving results table:
result.to_csv(r'Ma2Html_blitz.csv')

# plot the enrichment results and save to pdf
classes = ("CELL_FACTOR","CELL_GROWTH","CELL_WALL","FERMENTATION","MEMBRANE_TRANSPORT","METABOLISM_AMINO_ACID","METABOLISM_CARBON","METABOLISM_COFACTOR","METABOLISM_LIPIDS","METABOLISM_NUCLEIC_ACID","MOBILE_ELEMENT","MOBILITY","OPERONS_AUTOMATIC","REGULATIONS","RESPIRATION_ANAEROBIC","SECRETION","SPORULATION","STRESS","TRANSLATION","VIRULENCE_FACTORS")
nclasses = len(classes)

for i in range(nclasses):
   fig = blitz.plot.running_sum(signature, classes[i], library, result=result, compact=False)
   fig.savefig("running_sum_"+classes[i]+".svg", bbox_inches='tight')
   fig.savefig("running_sum_"+classes[i]+".tiff", bbox_inches='tight')
   plt.close(fig)

#fig_table = blitz.plot.top_table(signature, library, result, n=nclasses)
fig_table = blitz.plot.top_table(signature, library, result, n=9)
fig_table.savefig("top_table.svg", bbox_inches='tight')
fig_table.savefig("top_table.png", bbox_inches='tight')
fig_table.savefig("top_table.tiff", bbox_inches='tight')
