
import matplotlib.pyplot as pyp
import pandas as pd
from PIL import Image
import xlsxwriter
import pickle
import copy
from scipy.sparse import bsr_matrix
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from plotly.graph_objs import graph_objs
import plotly.figure_factory as ff
import plotly
from plotly.graph_objs import Scatter, Layout
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist





from collections import OrderedDict

from plotly import exceptions, optional_imports
from plotly.graph_objs import graph_objs

# Optional imports, may be None for users that only use our core functionality.
np = optional_imports.get_module('numpy')
scp = optional_imports.get_module('scipy')
sch = optional_imports.get_module('scipy.cluster.hierarchy')
scs = optional_imports.get_module('scipy.spatial')


def create_dendro(X, orientation="bottom", labels=None,
                      colorscale=None, distfun=lambda x: scs.distance.pdist(x),
                      linkagefun=lambda x: sch.linkage(x, method='complete'),
                      hovertext=None,color_threshold=None):
    """
    BETA function that returns a dendrogram Plotly figure object.
    :param (ndarray) X: Matrix of observations as array of arrays
    :param (str) orientation: 'top', 'right', 'bottom', or 'left'
    :param (list) labels: List of axis category labels(observation labels)
    :param (list) colorscale: Optional colorscale for dendrogram tree
    :param (function) distfun: Function to compute the pairwise distance from
                               the observations
    :param (function) linkagefun: Function to compute the linkage matrix from
                                  the pairwise distances
    :param (list[list]) hovertext: List of hovertext for constituent traces of dendrogram
        clusters
    Example 1: Simple bottom oriented dendrogram
    ```
    import plotly.plotly as py
    from plotly.figure_factory import create_dendrogram
    import numpy as np
    X = np.random.rand(10,10)
    dendro = create_dendrogram(X)
    plot_url = py.plot(dendro, filename='simple-dendrogram')
    ```
    Example 2: Dendrogram to put on the left of the heatmap
    ```
    import plotly.plotly as py
    from plotly.figure_factory import create_dendrogram
    import numpy as np
    X = np.random.rand(5,5)
    names = ['Jack', 'Oxana', 'John', 'Chelsea', 'Mark']
    dendro = create_dendrogram(X, orientation='right', labels=names)
    dendro['layout'].update({'width':700, 'height':500})
    py.iplot(dendro, filename='vertical-dendrogram')
    ```
    Example 3: Dendrogram with Pandas
    ```
    import plotly.plotly as py
    from plotly.figure_factory import create_dendrogram
    import numpy as np
    import pandas as pd
    Index= ['A','B','C','D','E','F','G','H','I','J']
    df = pd.DataFrame(abs(np.random.randn(10, 10)), index=Index)
    fig = create_dendrogram(df, labels=Index)
    url = py.plot(fig, filename='pandas-dendrogram')
    ```
    """
    if not scp or not scs or not sch:
        raise ImportError("FigureFactory.create_dendrogram requires scipy, \
                            scipy.spatial and scipy.hierarchy")

    s = X.shape
    if len(s) != 2:
        exceptions.PlotlyError("X should be 2-dimensional array.")

    if distfun is None:
        distfun = scs.distance.pdist

    dendrogram = _Dendrogram(X, orientation, labels, colorscale,
                             distfun=distfun, linkagefun=linkagefun,
                             hovertext=hovertext, color_threshold=color_threshold)

    return graph_objs.Figure(data=dendrogram.data, layout=dendrogram.layout)


class _Dendrogram(object):
    """Refer to FigureFactory.create_dendrogram() for docstring."""

    def __init__(self, X, orientation='bottom', labels=None, colorscale=None,
                 width="100%", height="100%", xaxis='xaxis', yaxis='yaxis',
                 distfun=lambda x: scs.distance.pdist(x),
                 linkagefun=lambda x: sch.linkage(x, method='complete'), #testing!
                 hovertext=None,color_threshold=None):
        self.orientation = orientation
        self.labels = labels
        self.xaxis = xaxis
        self.yaxis = yaxis
        self.data = []
        self.leaves = []
        self.sign = {self.xaxis: 1, self.yaxis: 1}
        self.layout = {self.xaxis: {}, self.yaxis: {}}

        if self.orientation in ['left', 'bottom']:
            self.sign[self.xaxis] = 1
        else:
            self.sign[self.xaxis] = -1

        if self.orientation in ['right', 'bottom']:
            self.sign[self.yaxis] = 1
        else:
            self.sign[self.yaxis] = -1

        if distfun is None:
            distfun = scs.distance.pdist

        (dd_traces, xvals, yvals,
            ordered_labels, leaves) = self.get_dendrogram_traces(X, colorscale,
                                                                 distfun,
                                                                 linkagefun, 
                                                                 hovertext,color_threshold)

        self.labels = ordered_labels
        self.leaves = leaves
        yvals_flat = yvals.flatten()
        xvals_flat = xvals.flatten()

        self.zero_vals = []

        for i in range(len(yvals_flat)):
            if yvals_flat[i] == 0.0 and xvals_flat[i] not in self.zero_vals:
                self.zero_vals.append(xvals_flat[i])

        self.zero_vals.sort()

        self.layout = self.set_figure_layout(width, height)
        self.data = graph_objs.Data(dd_traces)

    def get_color_dict(self, colorscale):
        """
        Returns colorscale used for dendrogram tree clusters.
        :param (list) colorscale: Colors to use for the plot in rgb format.
        :rtype (dict): A dict of default colors mapped to the user colorscale.
        """

        # These are the color codes returned for dendrograms
        # We're replacing them with nicer colors
        d = {'r': 'red',
             'g': 'green',
             'b': 'blue',
             'c': 'cyan',
             'm': 'magenta',
             'y': 'yellow',
             'k': 'black',
             'w': 'white'}
        default_colors = OrderedDict(sorted(d.items(), key=lambda t: t[0]))

        if colorscale is None:
            colorscale = [
                'rgb(0,116,217)',  # blue
                'rgb(35,205,205)',  # cyan
                'rgb(61,153,112)',  # green
                'rgb(40,35,35)',  # black
                'rgb(133,20,75)',  # magenta
                'rgb(255,65,54)',  # red
                'rgb(255,255,255)',  # white
                'rgb(255,220,0)']  # yellow

        for i in range(len(default_colors.keys())):
            k = list(default_colors.keys())[i]  # PY3 won't index keys
            if i < len(colorscale):
                default_colors[k] = colorscale[i]

        return default_colors

    def set_axis_layout(self, axis_key):
        """
        Sets and returns default axis object for dendrogram figure.
        :param (str) axis_key: E.g., 'xaxis', 'xaxis1', 'yaxis', yaxis1', etc.
        :rtype (dict): An axis_key dictionary with set parameters.
        """
        axis_defaults = {
                'type': 'linear',
                'ticks': 'outside',
                'mirror': 'allticks',
                'rangemode': 'tozero',
                'showticklabels': True,
                'zeroline': False,
                'showgrid': False,
                'showline': True,
            }

        if len(self.labels) != 0:
            axis_key_labels = self.xaxis
            if self.orientation in ['left', 'right']:
                axis_key_labels = self.yaxis
            if axis_key_labels not in self.layout:
                self.layout[axis_key_labels] = {}
            self.layout[axis_key_labels]['tickvals'] = \
                [zv*self.sign[axis_key] for zv in self.zero_vals]
            self.layout[axis_key_labels]['ticktext'] = self.labels
            self.layout[axis_key_labels]['tickmode'] = 'array'

        self.layout[axis_key].update(axis_defaults)

        return self.layout[axis_key]

    def set_figure_layout(self, width, height):
        """
        Sets and returns default layout object for dendrogram figure.
        """
        self.layout.update({
            'showlegend': False,
            'autosize': False,
            'hovermode': 'closest',
            'width': width,
            'height': height
        })

        self.set_axis_layout(self.xaxis)
        self.set_axis_layout(self.yaxis)

        return self.layout

    def get_dendrogram_traces(self, X, colorscale, distfun, linkagefun, hovertext,color_threshold):
        """
        Calculates all the elements needed for plotting a dendrogram.
        :param (ndarray) X: Matrix of observations as array of arrays
        :param (list) colorscale: Color scale for dendrogram tree clusters
        :param (function) distfun: Function to compute the pairwise distance
                                   from the observations
        :param (function) linkagefun: Function to compute the linkage matrix
                                      from the pairwise distances
        :param (list) hovertext: List of hovertext for constituent traces of dendrogram
        :rtype (tuple): Contains all the traces in the following order:
            (a) trace_list: List of Plotly trace objects for dendrogram tree
            (b) icoord: All X points of the dendrogram tree as array of arrays
                with length 4
            (c) dcoord: All Y points of the dendrogram tree as array of arrays
                with length 4
            (d) ordered_labels: leaf labels in the order they are going to
                appear on the plot
            (e) P['leaves']: left-to-right traversal of the leaves
        """
        d = distfun(X)
        Z = linkagefun(d)
        P = sch.dendrogram(Z, orientation=self.orientation,
                           labels=self.labels, no_plot=True,color_threshold=color_threshold)

        icoord = scp.array(P['icoord'])
        dcoord = scp.array(P['dcoord'])
        ordered_labels = scp.array(P['ivl'])
        color_list = scp.array(P['color_list'])
        colors = self.get_color_dict(colorscale)

        trace_list = []

        for i in range(len(icoord)):
            # xs and ys are arrays of 4 points that make up the '∩' shapes
            # of the dendrogram tree
            if self.orientation in ['top', 'bottom']:
                xs = icoord[i]
            else:
                xs = dcoord[i]

            if self.orientation in ['top', 'bottom']:
                ys = dcoord[i]
            else:
                ys = icoord[i]
            color_key = color_list[i]
            hovertext_label = None
            if hovertext:
                hovertext_label = hovertext[i]
            trace = graph_objs.Scatter(
                x=np.multiply(self.sign[self.xaxis], xs),
                y=np.multiply(self.sign[self.yaxis], ys),
                mode='lines',
                marker=graph_objs.Marker(color=colors[color_key]),
                text=hovertext_label,
                hoverinfo='text'
            )

            try:
                x_index = int(self.xaxis[-1])
            except ValueError:
                x_index = ''

            try:
                y_index = int(self.yaxis[-1])
            except ValueError:
                y_index = ''

            trace['xaxis'] = 'x' + x_index
            trace['yaxis'] = 'y' + y_index

            trace_list.append(trace)

        return trace_list, icoord, dcoord, ordered_labels, P['leaves']



class Dendrogram:
    def __init__(self,file,min_d=0.8E-3,**args):
        '''
        Initializes BioDendro - package used to cluster and plot dendrograms
    
        '''
        
        #Converts **args to self.variable
        tmp=copy.deepcopy(locals()) #Don't copy list to another list!
        for each in args:
            if not "__" in each:
                exec ("self."+each+"="+str(args[each]))
        if 'sheetname' in args:
            myfile=pd.read_excel(file,sheetname=int(self.sheetname)) #pandas df that reads the file
        else:
            myfile=pd.read_excel(file)

        self.CUTOFF=5.0
        self.myfile=myfile
        self.GL=False
        self.CLUSTERIZE=False
    
    def clusterize(self):
        '''
        Clusters based on a metric passed as input. Defaults are based on data in the
        test csv file
        Input: csv
        Outputs: 
        '''
        #Default parameters - column name, sheet name, minimum distance to cluster
        colname='mz' #Look for a particular column to cluster
        sheetname=1 #Which sheet in the notebook?
        min_d=0.8E-3 #Cluster cutoff used
        myfile=self.myfile
        clusters=[]
        clusters.append(0)
        clusters=list(myfile[myfile[colname].diff()>=min_d].index) #Label based on cluster cutoff
        myfile['labels']=0
        cnt=1
        col_names=[]
        new_col=np.zeros(len(myfile),dtype=np.int)
        clusters.insert(0,0)
        #The routine below generates 
        for each in clusters:
            if cnt==1:
                for tmp_cnt in range(cnt-1,cnt+1):
                    new_col[0:clusters[tmp_cnt]+1]=cnt
                    if clusters[tmp_cnt]!=0:
                        tmp=myfile['mz'][0:clusters[tmp_cnt]+1]
                    else:
                        tmp=myfile['mz'][0]
                    col_names.append(str("%.4f"%np.around(np.mean(tmp),4))+
                                      '_'+str("%.4f"%np.around(np.min(tmp),4))+
                                      '_'+str("%.4f"%np.around(np.max(tmp),4)))
            elif cnt==len(clusters):
                new_col[clusters[cnt-1]:]=cnt
                tmp=myfile['mz'][clusters[cnt-1]:]
                col_names.append(str("%.4f"%np.around(np.mean(tmp),4))+
                                 '_'+str("%.4f"%np.around(np.min(tmp),4))+
                                 '_'+str("%.4f"%np.around(np.max(tmp),4)))
            else:
                new_col[clusters[max(cnt-1,0)]:clusters[cnt]]=cnt
                tmp=myfile['mz'][clusters[max(cnt-1,0)]:clusters[cnt]]
                col_names.append(str("%.4f"%np.around(np.mean(tmp),4))+
                                 '_'+str("%.4f"%np.around(np.min(tmp),4))+
                                 '_'+str("%.4f"%np.around(np.max(tmp),4)))

            cnt+=1
        #Check for consistency of labels #0, and #1 as diff condition requires that the first 
        #label is sacrificed for binning.
        #clusters->starting points of clusters, new_col->cluster id
        #Compensate for the diff
        myfile['labels']=list(new_col)
        self.myfile=myfile
        self.clusters=clusters
        self.new_col=new_col
        self.cnt=cnt
        self.col_names=col_names
        self.indices=myfile['sample']
	self.full_df=""  #DataFrame with bin ids as columns and records as rows
        

        #Flags set to see to check and see what functions are executed
        self.GL=False #Generate Linkage
        self.GO=False #Generate output
        self.DENDRO=False #Check if dendro was created
        self.CLUSTERIZE=True
        self.LINKAGE=False
        #return myfile,clusters,new_col,cnt,col_names

    def pop_filled_matrices(self,matrix):
        '''Returns the column ids that are similar'''
        tmp=matrix.any(axis=0)
        cnt=0
        filled_indices=[]
        for each in tmp:
            if (each and cnt<matrix.shape[1]-2):
                filled_indices.append(cnt)
            cnt+=1
        return filled_indices
    

    def get_indices(self,tf,indices):
            '''What are the index ids??'''
            cnt=0
            mylist=[]
            total_cnt=0
            for each in tf:
                    if each :
                            mylist.append(indices[cnt])
                            total_cnt+=1
                    cnt+=1
            return mylist
        

    def plot_bins(self,inp,filename):
        vals=np.sum(inp,axis=0)/len(inp)
        pyp.clf()
        pyp.bar(range(len(vals)),vals)
        pyp.savefig(filename)

    def generate_linkage(self):
        A=[]
        cluster=dict()
        cnt=0
        max_val=np.max(self.myfile['labels'])+1
        for each in self.myfile['sample']:
            tmpval=each
            loc=self.myfile['labels'].iloc[cnt]
            if tmpval in cluster:
                #if loc>10:
                cluster[tmpval][loc]=True

            else:
                cluster[tmpval]=np.zeros(max_val,dtype=np.bool)
                cluster[tmpval][loc]=False
            A.append(tmpval)
            cnt+=1
        self.cluster=cluster
        Full_matrix=np.zeros([len(cluster),max_val],dtype=np.bool)
        cnt=0
        mymax=0
        self.labels=[]
        for key in cluster:
            Full_matrix[cnt,:]=cluster[key]
            self.labels.append(key)
            npsum=np.sum(cluster[key])
            if npsum>mymax:
                mymax=npsum
            cnt+=1
        self.full=Full_matrix
        if not self.LINKAGE:
            self.Dis=pdist(self.full)
            self.Z=linkage(self.Dis, method='complete')
            self.LINKAGE=True
        self.mycluster=fcluster(self.Z,self.CUTOFF,criterion='distance')
        self.GL=True


    def write_full_matrix(self,filename="full_matrix.xlsx"):
	if self.GL==True:
        	df=pd.DataFrame(self.full*1)
        	df.columns=self.col_names
        	df['index']=self.labels
        	df.set_index(['index'])
        	excel_writer=pd.ExcelWriter(filename)
        	df.to_excel(excel_writer,'Sheet 1')
		self.full_df=df

        
    def generate_out(self):
        cnt=0
        colnames=pd.DataFrame(self.col_names)
        for each in range(np.max(self.mycluster)+1):
        #for each in range(1,np.max(self.mycluster)-1):
            tmp_indices=self.get_indices(self.mycluster==each,self.indices)
            with open("results/cluster_"+str(len(tmp_indices))+"_"
                     +str(cnt)+".txt","a") as text_file:
                tmp_mat=pd.DataFrame(self.full[self.mycluster==each])
                myiloc=self.pop_filled_matrices(self.full[self.mycluster==each])
                tmp_indices=self.get_indices(self.mycluster==each,self.labels)
                mytmp=pd.DataFrame(tmp_mat*1)
                mytmp.index=tmp_indices
                mytmp=mytmp.iloc[:,myiloc]
                mytmp.columns=colnames.iloc[myiloc]
                cluster_label="Cluster : "+str(cnt+1)+" Length : "+str(len(tmp_indices))+"\n"
                text_file.write(cluster_label)
                mytmp.to_csv(text_file,sep="\t")
                self.plot_bins(mytmp,"results/Cluster_"+str(cnt)+".png")
                cnt+=1
        print (cnt)
        self.GO=True


    def visualize(self,cutoff=5.0,x=900,y=600):
    #Author: Paula
        if not self.CUTOFF==cutoff:
           self.GL=False
        if not self.CLUSTERIZE:
           self.clusterize()
        if not self.GL:
           self.generate_linkage()
        if not self.GO:
           self.generate_out()
        if self.GL==True: #Check whether linkage and outputs are generated
           c,coph_dist=cophenet(self.Z,pdist(self.full))
           rnames=list(self.indices)
           self.dendro=create_dendro(self.full,labels=rnames,color_threshold=cutoff)
           self.dendro['layout'].update({'width':x, 'height':y, 'title':'Python plotly dendro', 'xaxis':{'title':'sample'}, 'yaxis':{'title':'distance'}, 'hovermode':'closest'})
           self.dendro['data'].update({'hoverinfo':'all'})
           plotly.offline.plot(self.dendro,filename='simple_dendrogram.html')
        return self.dendro
