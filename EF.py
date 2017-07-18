'''Enrichment Factor Calculation'''
'''
Calculates Enrichment Factor Values for the input maps given
'''

#Modules
import numpy as np

#Function
def EF(luc,d,cdl,cd,N,omap,amap,mask,row,column):
    #Transition analysis dictionary (used to store presence count) 
    #Determine the composition of the neighbourhood of cells that transitioned
    #between the time steps from the data
    transition_dictionary={}
    for p in range(0,luc):
        for q in range(0,luc):
            for c in range(0,d):
                key="c-"+str(p)+"|n-"+str(q)+"|d-"+str(c)
                value=[0]*(1+N[c])
                transition_dictionary[key]=value
    #Intialise a numpy array to store the sum of enrichments for neighbourhood analysis
    EF_sum=np.zeros(shape=(d,luc,luc))
    #Initalise a list to store the count of transitions from the data and the number of luc
    apl_count=[0]*luc
    luc_count=[0]*luc
    #Extract a count of the presence of different land-use classes in neighbourhooods
    for i in range(0,row):
        for j in range(0,column):
            #Skip if masked out of the region map
            if mask[i,j]<1:                     
                next
            elif omap[i,j]>luc-1:
                next
            elif amap[i,j]>luc-1:
                next
            else:
                #Add count
                luc_count[omap[i,j]]=luc_count[omap[i,j]]+1
                #Peform analysis if cells have transitioned
                if omap[i,j]==amap[i,j]:
                    next
                else:
                    #Initialise an array to store the count of different land-use classes in the neighbourhood of the cell of interest      
                    float_store_count=np.zeros(shape=(luc,cdl))
                    neighbourhood_size=[0]*d            
                    for c in range(0,d+1):
                        #Analyse the presence at the location of interest
                        if c==0:
                            k=i
                            l=j
                            if mask[k,l]==0:
                                next
                            else:
                                float_store_count[omap[k,l],0]=float_store_count[omap[k,l],0]+1    
                        #Iterate around the neighbourhood of the cell of interest 
                        else:
                            k=i-c
                            if k<0:
                                next
                            else:
                                for l in range(j-c,j+c+1):
                                    if l<0:
                                        next
                                    elif l>=column:
                                        next
                                    elif mask[k,l]==0:
                                        next
                                    else:
                                        idx=((k-i)**2+(l-j)**2)**0.5
                                        for e in range(0,cdl):
                                            if idx==cd[e]:
                                                float_store_count[omap[k,l],e]=float_store_count[omap[k,l],e]+1
                            k=i+c
                            if k>=row:
                                next
                            else:
                                for l in range(j-c, j+c+1):
                                    if l<0:
                                        next
                                    elif l>=column:
                                        next
                                    elif mask[k,l]==0:
                                        next
                                    else:
                                        idx=((k-i)**2+(l-j)**2)**0.5
                                        for e in range(0,cdl):
                                            if idx==cd[e]:
                                                float_store_count[omap[k,l],e]=float_store_count[omap[k,l],e]+1
                            l=j-c
                            if l<0:
                                next
                            else:
                                for k in range(i-c+1,i+c):
                                    if k<0:
                                        next
                                    elif k>=row:
                                        next
                                    elif mask[k,l]==0:
                                        next
                                    else:
                                        idx=((k-i)**2+(l-j)**2)**0.5
                                        for e in range(0,cdl):
                                            if idx==cd[e]:
                                                float_store_count[omap[k,l],e]=float_store_count[omap[k,l],e]+1
                            l=j+c
                            if l>=column:
                                next
                            else:
                                for k in range(i-c+1,i+c):
                                    if k<0:
                                        next
                                    elif k>=row:
                                        next
                                    elif mask[k,l]==0:
                                        next
                                    else:
                                        idx=((k-i)**2+(l-j)**2)**0.5
                                        for e in range(0,cdl):
                                            if idx==cd[e]:
                                                float_store_count[omap[k,l],e]=float_store_count[omap[k,l],e]+1           
                    #Initialise an array to aggregate the different discrete distances into rings 1, 2, 3
                    #Rings aggregate based on proximity. 0.5->1<-1.5, 1.5->2<-2.5 etc.            
                    float_store_aggregated=np.zeros(shape=(luc,d))
                    neighbourhood_size=[0]*d     
                    for p in range(0,luc):
                        for c in range(0,cdl):
                            if float_store_count[p,c]>0:
                                for e in range(0,d):
                                    x=cd[c]
                                    if (e-0.5)<x<(e+0.5):
                                        float_store_aggregated[p,e]=float_store_aggregated[p,e]+float_store_count[p,c]
                                        neighbourhood_size[e]=neighbourhood_size[e]+float_store_count[p,c]
                    #Convert float_store_aggregated to proportion values for processing (Range 0-1, neighbourhood size dependent)
                    fsa_proportion=np.zeros(shape=(luc,d))
                    for p in range(0,luc):
                        for c in range(0,d):
                            if float_store_aggregated[p,c]==0:
                                next
                            else:
                                fsa_proportion[p,c]=float_store_aggregated[p,c]/neighbourhood_size[c]
                    #Store values into transition dictionary bins (widths .5/total number) 
                    central=amap[i,j]
                    apl_count[central]=apl_count[central]+1
                    for p in range(0,luc):
                        for c in range(0,d):
                            EF_sum[c,central,p]=EF_sum[c,central,p]+fsa_proportion[p,c]
                            key="c-"+str(central)+"|n-"+str(p)+"|d-"+str(c)
                            for a in range(0,N[c]+1):
                                lower_bound=float(a)/float(N[c])-(0.5/N[c])
                                upper_bound=float(a)/float(N[c])+(0.5/N[c])
                                if lower_bound<fsa_proportion[p,c]<upper_bound:
                                    transition_dictionary[key][a]=transition_dictionary[key][a]+1
    '''Process enrichment factor calcualtion'''
    EF=np.zeros(shape=(d,luc,luc))
    for p in range(0,luc):
        for q in range(0,luc):
            for c in range(0,d):
                if apl_count[p]>0 and luc_count[q]>0:
                    x=EF_sum[c,p,q]/apl_count[p]
                    y=float(luc_count[q])/sum(luc_count[:])
                    EF[c,p,q]=x/y
    return EF
                            