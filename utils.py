import pandas as pd
from pandas.api.types import CategoricalDtype
import numpy as np
from math import *
from collections import Counter

# wimport scanpy as sc
# import scvelo as scv

import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

def run(cmd):
    proc = subprocess.Popen([cmd], shell=True,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE,
    )
    stdout, stderr = proc.communicate()
    return proc.returncode, stdout, stderr

def getarg(argv, index, default=None):
    if len(argv) > index:
        ret = argv[index]
    else:
        ret = default
    return ret

def read_gmt(fpath):
    gdict = {}
    
    with open(fpath, 'r') as f:
        lines = f.readlines()
        for line in lines:
            items = line.strip().rstrip().split('\t')
            gset = items[0].lower()
            genes = items[2:]
            gdict[gset] = genes            
    return gdict

def write_gmt(gdict, fpath):  
    with open(fpath, 'w') as f:
        for key, val in gdict.items():
            val_string = '\t'.join(val)
            f.write('%s\thttp\t%s\n'%(key, val_string))            
    print('finished writing %d gene sets to %s'%(len(gdict), fpath))

    
def plot_by_single_gene2(genes, df_p,ncols, nrows,group='class_2c'):
    
    def join_with_linebreaks(words):
        modified_words = [word if i%4 != 3 else word + '\n' for i, word in enumerate(words)]
        return ' '.join(modified_words)

    from statannotations.Annotator import Annotator
    
    class_order = ['Non-proliferative', 'Proliferative']
    pairs = [(class_order[i],class_order[j]) \
             for i in range(len(class_order)) for j in range(len(class_order)) if i<j ]

    sns.set(font_scale=1.7)
    sns.set_style('white')

    fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 6*nrows))
    plt.subplots_adjust(wspace=3, hspace=5)
    
    my_palette = {'Non-proliferative': 'lightblue', 'Proliferative': 'orange'}

    for i, gene in enumerate(genes):
        ax = axes.flatten()[i]; gene_data = df_p.loc[:,[gene,group]]

        sns.boxplot(data=gene_data, x=group,y=gene, ax=ax, palette=my_palette, order=class_order)
        sns.swarmplot(data=gene_data, x=group,y=gene, ax=ax, color='k', alpha=0.7, order=class_order)

        title = join_with_linebreaks(gene.split('_'))

        ax.set_title('ssGSEA Score')
        
        
        ax.set_xticklabels(['Hormone/\nImmune-related', 'Proliferative'], rotation=0)
        ax.set_xlabel('')
        ax.set_ylabel(title)
        
        annotator = Annotator(ax, pairs, data=df_p, x=group, y=gene, 
                          comparisons_correction='fdr_bh')
        annotator.configure(test='Mann-Whitney', text_format='star', loc='inside')
        annotator.apply_and_annotate()

    for i in range(nrows*ncols):
        ax = axes.flatten()[i]
        if i >= len(genes):
            ax.axis('off')
            
    plt.tight_layout()
    plt.show()

