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

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


#data refactoring with pandas
def get_row_data():
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
	return data

row_data = get_row_data()
data_full = pandas.DataFrame(data = row_data, columns = ['Door1','Door2','FaucetHandle1','FaucetHandle2','SinkFloor1','SinkFloor2','Soap Dispenser','Stall','ToiletFloor1','ToiletFloor2','ToiletFlushHandle1','ToiletFlushHandle2','ToiletSeat1','ToiletSeat2'])

#Contruction du dataset pour mieux visualiser des corrélations
def get_dataset(data_set, section, threshold):
	AMB = pandas.pivot_table(data_set, columns = section, aggfunc='sum')
	vec = AMB.apply(lambda col:col.sum()>=threshold, axis=0)

	AMB=AMB.loc[:,vec]
	data = pandas.DataFrame.transpose(AMB)
    
	return data


#creation du jeu de données avec embranchement threshold de 10
Section = 'Embranchement'
data = get_dataset(row_data, Section, 10)



# figure clustering
def get_hac(data, label):
	Z = hac.linkage(data, method='complete', metric='correlation', optimal_ordering=False)

    # Plot dendogram
	plt.figure(figsize=(18, 10))
	plt.title('Hierarchical Clustering Dendrogram')
	plt.xlabel('Sample name')
	plt.ylabel('Distance')
	hac.dendrogram(
           Z,
           leaf_rotation=45.,  # rotates the x axis labels
           leaf_font_size=12.,  # font size for the x axis labels
           labels = label
	)
	plt.tight_layout()
	plt.savefig('fig/dendogram.png')


get_hac(data_full.T,data_full.columns)

# Figure : PCA
def get_PCA(data,Section):
	data_PCA = data.values
	data_PCA = StandardScaler().fit_transform(data_PCA)
	pca = PCA(n_components=2)
	principalComponents = pca.fit_transform(data_PCA)


	principalDf = pandas.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2'])
	tagretDf = pandas.DataFrame(data = data.columns, columns = ['target'])
	finalDf = pandas.concat([principalDf, tagretDf[['target']]], axis = 1)

	fig = plt.figure(figsize = (8,8))
	ax = fig.add_subplot(1,1,1) 
	ax.set_xlabel('Principal Component 1 : {} %'.format(round(pca.explained_variance_ratio_[0]*100,2)), 
                  fontsize = 15, labelpad = 15)
	ax.set_ylabel('Principal Component 2 : {} %'.format(round(pca.explained_variance_ratio_[1]*100,2)), 
                  fontsize = 15, labelpad = 15)
	ax.set_title('2 Component PCA : {}'.format(Section), fontsize = 20, pad = 20)
    
	for target in data.columns:
		indicesToKeep = finalDf['target'] == target
		ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1'], 
                   finalDf.loc[indicesToKeep, 'principal component 2'], 
                   c = 'black', s = 30 )

		plt.annotate(target, # this is the text
                 (finalDf[indicesToKeep]['principal component 1'],
                  finalDf[indicesToKeep]['principal component 2']), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center', size=10) # horizontal alignment can be left, right or center
		plt.tight_layout()
		plt.savefig('fig/pca.png')
		
get_PCA(data, Section)


#Figure heatmap
def get_heatmap(data, section):
    # get correlations
	fig, ax = plt.subplots(figsize=(12, 10))

    #corelation
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
               linewidths=0.3, vmin=0, vmax=1, cbar=False,
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

get_heatmap(data, Section)


#Figure 2.1: clustering hierarchique
def get_clustering(data_set, section):
    # color map
	cmap = sns.diverging_palette(0, 230, 90, 60, as_cmap=True)
    
    #plot clustermap
	sns.clustermap(data_set, cmap="Blues", vmin= 0, vmax=2800,metric='euclidean', #euclidean
                   figsize=(10, 10),dendrogram_ratio=0.3,colors_ratio=0.003,cbar_pos=(0.02, 0.8, 0.05, 0.18),
                   linewidth=0.3, cbar_kws={"shrink": .8})
    
    # axis labels
	plt.xlabel('Sample')
	plt.ylabel(section)

    # title
	plt.title('Expression level : {} / Sample'.format(section), loc='left')
	plt.tight_layout()
	plt.savefig('fig/hicluster.png')

get_clustering(data, Section)


#Figure Sunburst
#modéliser la présence des espèces par leur expression dans chacune des conditions
data = get_row_data()

sample='FaucetHandle'

def make_figure(id_sample):
    # make sunburst
	sb1 = px.sunburst(data, path=['Domaine','Embranchement','Classe','Ordre','Famille','Genre'], 
                      values=round((data[str(id_sample)+"1"]/sum(data[str(id_sample)+"1"]))*100,2),
                      color='Embranchement', branchvalues='total',template='ggplot2',title='sunburst',
                      color_discrete_sequence=px.colors.qualitative.Pastel)._data

	sb2 = px.sunburst(data, path=['Domaine','Embranchement','Classe','Ordre','Famille','Genre'], 
                      values=round((data[str(id_sample)+"2"]/sum(data[str(id_sample)+"2"]))*100,2),
                      color='Embranchement', branchvalues='total',template='ggplot2',title='sunburst',
                      color_discrete_sequence=px.colors.qualitative.Pastel)._data

    # traces with separate domains to form a subplot
	trace1 = go.Sunburst(ids=sb1[0]['ids'],
                         labels=sb1[0]['labels'],
                         values=sb1[0]['values'],
                         parents=sb1[0]['parents'],
                         domain={'x': [0.0, 0.5], 'y': [0, 1]},
                         branchvalues ="total")

	trace2 = go.Sunburst(ids=sb2[0]['ids'],
                         labels=sb2[0]['labels'],
                         values=sb2[0]['values'],
                         parents=sb2[0]['parents'],
                         domain={'x': [0.5, 1.0], 'y': [0, 1]},
                         branchvalues ="total")

    # layout and figure production
	layout = go.Layout(height = 500,
                       width = 950,
                       autosize = False,
                       title = 'Side by side Sunburst diagrams by kind of paired sample : Sample 1 vs Sample 2')

	fig = go.Figure(data = [trace1, trace2], layout = layout)
	fig.update_traces(textinfo='label+percent entry')
	fig.write_image('fig/sunbrust.png')

make_figure(sample)

