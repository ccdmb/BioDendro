import numpy as np
import matplotlib.pyplot as pyp
import pandas as pd
from PIL import Image
import xlsxwriter
import pickle
import copy
from scipy.sparse import bsr_matrix
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster



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

        self.myfile=myfile
    
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
                new_col[0:clusters[cnt]]=cnt
                tmp=myfile['mz'][0:clusters[cnt]]
                col_names.append(str("%.4f"%np.around(np.mean(tmp),4))+
                                 '_'+str("%.4f"%np.around(np.min(tmp),4))+
                                 '_'+str("%.4f"%np.around(np.max(tmp),4)))
                #col_names.append(str(np.mean(myfile['mz']))+'_'+str(np.std(inp))+'_'+str(np.min(inp))+'_'+str(np.max(inp))clust)
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
        #clusters->starting points of clusters, new_col->cluster id
        #Compensate for the diff
        myfile['labels']=list(new_col)
        self.myfile=myfile
        self.clusters=clusters
        self.new_col=new_col
        self.cnt=cnt
        self.col_names=col_names
        
        return myfile,clusters,new_col,cnt,col_names

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
                    if each:
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
        for key in cluster:
            Full_matrix[cnt,:]=cluster[key]
            npsum=np.sum(cluster[key])
            if npsum>mymax:
                mymax=npsum
            cnt+=1
        self.full=Full_matrix
        self.Z=linkage(self.full)
        self.mycluster=fcluster(self.Z,5.0,criterion='distance')
        
    def generate_out(self):
        cnt=0
        colnames=pd.DataFrame(self.col_names)
        for each in range(1,np.max(self.mycluster)-1):
            tmp_indices=self.get_indices(self.mycluster==each,self.col_names)
            with open("results/cluster_"+str(len(tmp_indices))+"_"
                     +str(cnt)+".txt","a") as text_file:
                tmp_mat=pd.DataFrame(self.full[self.mycluster==each])
                myiloc=self.pop_filled_matrices(self.full[self.mycluster==each])
                tmp_indices=self.get_indices(self.mycluster==each,self.col_names)
                mytmp=pd.DataFrame(tmp_mat*1)
                mytmp.index=tmp_indices
                mytmp=mytmp.iloc[:,myiloc]
                mytmp.columns=colnames.iloc[myiloc]
                cluster_label="Cluster : "+str(cnt+1)+" Length : "+str(len(tmp_indices))+"\n"
                text_file.write(cluster_label)
                mytmp.to_csv(text_file,sep="\t")
                self.plot_bins(mytmp,"results/Cluster_"+str(cnt)+".png")
                cnt+=1
