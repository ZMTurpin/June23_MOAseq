#import necessary packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob

#read in bedpe files

files = glob.glob('*bedpe')
files.sort()

names = []
for a in files:
    b=a.split("_")
    names.append(b[0]+"_"+b[1]+"_"+b[2])

beds = []

for a in files:
    beds.append(pd.read_csv(a, sep='\t', header=None, index_col=None))
#compute size column
for b in beds:
    b['size'] = b[2] - b[1]
#plot a histogram of 'size' for each dataframe

fig,axs = plt.subplots(nrows=len(beds), ncols=1, sharex=True, sharey=True, figsize=(4,8))

for c, ax in enumerate(axs.flatten()):
    ax.hist(beds[c]['size'], bins=50)
    ax.set_title(names[c])
    ax.set_xlabel('Aligned Fragment Size (bp)', fontsize=12)
    ax.set_ylabel('Frequency')
    ax.axvline(x=beds[c]['size'].mean(), ymin=0, ymax=1, linestyle="--", color='red')

plt.tight_layout()
plt.savefig('SeedlingFlooding_MOA_bedHist.pdf', format='pdf')
