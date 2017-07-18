'''                     CONSIDERED DISTANCES                     '''
''' Determines how many delineations of neighbourhoods is requried 
based on the maximum distance considered specificed '''

def considered_distances(d):
    considered_distances=[]
    for i in range (1,d):                                             # Create a list of distances that must be considered
        for j in range (0,d):                                         # based on square arrangments. This WONT include the
            x = (i**2+j**2)**0.5                                      # maximum distance, which is added below
            if x in considered_distances:
                next
            elif x>=d:
                next
            else:
                considered_distances.append(x)
    cdl=len(considered_distances)                                     # Find the number of considered distances
    considered_distances.append(d)
    considered_distances.sort()                                       # Sort values smallest to largest
    considered_distances=[0]+considered_distances
    cdl=len(considered_distances)                                     # Find the number of considered distances
    return considered_distances,cdl