class Dendrogram:
    def __init__(self,file,min_d=0.8E-3,**args):
        '''
        Initializes BioDendro - package used to cluster and plot dendrograms
    
        '''
        
        #Converts **args to self.variable
        tmp=copy.deepcopy(locals()) #Don't copy list to another list!
        for each in tmp:
            if not "__" in each:
                exec ("self."+each+"="+each)
        myfile=pd.read_excel(file,sheetname=1) #pandas df that reads the file
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
        self.myfile=myfile
        self.clusters=clusters
        self.new_col=new_col
        self.cnt=cnt
        self.col_names=col_names
        
        return myfile,clusters,new_col,cnt,col_names
