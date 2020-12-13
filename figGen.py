import matplotlib.pyplot as plt
import pandas
import seaborn as sns
import scipy
import numpy as np

import scipy.cluster.hierarchy as hac

import plotly
import plotly.offline as py
import plotly.graph_objs as go
import plotly.express as px

from jupyter_dash import JupyterDash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

level = 'Embranchement'

#data refactoring with pandas
data = pandas.read_table('1335.vdb.tab', sep = '\t', header = 0)
data = data.sort_values(by = ['OTUConTaxonomy'])
data = data.reset_index(drop=True)

list_taxon=['Domaine','Embranchement','Classe','Ordre',"Famille","Genre"]
for i in range(6):
    colname = list_taxon[i]
    data[colname]=data.apply(lambda row:row.OTUConTaxonomy.split(";")[i].split("(")[0],axis=1)
    
data['OTU'] = data.apply(lambda row:row.OTUNumber.replace('Otu',''), axis=1)
data = data.drop(columns=['OTUConTaxonomy','OTUNumber'])
data.columns=['Door1','Door2','FaucetHandle1','FaucetHandle2','SinkFloor1','SinkFloor2','Soap Dispenser','Stall','ToiletFloor1','ToiletFloor2','ToiletFlushHandle1','ToiletFlushHandle2','ToiletSeat1','ToiletSeat2','Domaine','Embranchement','Classe','Ordre','Famille','Genre','OTU']

#Figure 1 : barplot
#EMB = data.groupby(by=['Embranchement']).sum()

AMB = pandas.pivot_table(data,columns = level,aggfunc='sum')
vec = AMB.apply(lambda col:col.sum()>=10, axis=0)


AMB=AMB.loc[:,vec]

fig, ax = plt.subplots() 
ax = AMB.plot(kind="bar",stacked=True)
ax.legend(loc='lower left',fontsize=10, bbox_to_anchor=(1.1,0))
ax.set_ylabel("Aboundance [nbr de lectures]")
ax.set_title("Bar plot des {}".format(level))

plt.tight_layout()
plt.savefig('fig/barplot.png')

#Contruction du dataset pour les heatmaps
def get_dataset(data_set, section, threshold):
    AMB = pandas.pivot_table(data_set,columns = section,aggfunc='sum')
    vec = AMB.apply(lambda col:col.sum()>=threshold, axis=0)

    AMB=AMB.loc[:,vec]
    data_hm = pandas.DataFrame.transpose(AMB)
    return data_hm

data_hm = get_dataset(data, level, 10)

#Figure 2.1: clustering hierarchique
def get_clustering(data_set, section):
    #fig, ax = plt.subplots(figsize=(12, 10))
    
    # color map
	cmap = sns.diverging_palette(0, 230, 90, 60, as_cmap=True)
    
	sns.clustermap(data_set, cmap="Blues", vmin= 0, vmax=2800,metric='euclidean', #euclidean
                   figsize=(10, 10),dendrogram_ratio=0.3,colors_ratio=0.003,cbar_pos=(0.02, 0.8, 0.05, 0.18),
                   linewidth=0.3, cbar_kws={"shrink": .8})
    
    # axis labels
	plt.xlabel('Sample')
	plt.ylabel(section)

	# title
	plt.title('Expression level : {} / Sample'.format(section), loc='left')
	plt.tight_layout()
	plt.savefig('fig/cluster.png')

get_clustering(data_hm, level)

#Figure 2.2: heatmap
def get_heatmap(data, section):
    # get correlations
	fig, ax = plt.subplots(figsize=(12, 10))

    #corelation
    #print(data_hm)
	df_corr = data.corr()# irrelevant fields
	np.ones_like(df_corr, dtype=np.bool)

    # mask
	mask = np.triu(np.ones_like(df_corr, dtype=np.bool))

    # adjust mask and df
	mask = mask[1:, :-1]
	corr = df_corr.iloc[1:,:-1].copy()

    # color map
	cmap = sns.diverging_palette(0, 230, 90, 60, as_cmap=True)


    # plot heatmap
	sns.heatmap(corr, mask=mask, annot=True, fmt=".2f", cmap="Blues",
               linewidths=0.3, vmin=-1, vmax=1, cbar=False,
               cbar_kws={"shrink": .8}, square=True)

    # ticks
	yticks = [i for i in corr.index]
	xticks = [i for i in corr.columns]

	plt.yticks(plt.yticks()[0], labels=yticks, rotation=0)
	plt.xticks(plt.xticks()[0], labels=xticks)

    # title
	plt.title('Expression level : Sample / Sample for the section "{}"'.format(section), loc='center')
	plt.tight_layout()
	plt.savefig('fig/heatmap.png')

get_heatmap(data_hm, level)

