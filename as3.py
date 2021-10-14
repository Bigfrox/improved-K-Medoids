'''
Data Mining Assignment 3, Improved K-Medoids
2016253072
명수환(Myeong Suhwan)
'''
from datetime import datetime
import math

from numpy import short



'''
Improved k-Medoids Algorithm

Step1. Initial medoids selection:
        1. Calculate the distance between every pair of objects.
        2. Calculate the sum of distance for each object.
        3. Select k objects having the smallest sum of distance as initial medoids.
        4. Obtain the initial clusters by assigning each non-medoid to the nearest medoid.

Step2. Update medoids iteratively:
        1. For each cluster, calculate the sum of distance within the cluster for each object and select a new medoid having the smallest sum of distance.
        2. Obtain the updated clusters by assigning each mon-medoid to the nearest medoid.
        3. Repeat steps 2-1 and 2-2 until the clusters do not change.
'''

def getDataFromFile(filename):
    input_file = open(filename, 'r')
    gene_id = dict()
    
    line_num = 0
    for line in input_file:
        time_point_data = line.split()
        gene_id[line_num] = time_point_data
        
        line_num += 1
        
        
    return gene_id

def randClustering(gene_id):
    cluster = list() # k개의 cluster 생성, k : 10
    init_mem_num_in_cluster = int(len(gene_id) / k)
    line_num = 0

    for i in range(0,len(gene_id)-1,init_mem_num_in_cluster):
        same_cluster = list()
        for i in range(init_mem_num_in_cluster):
            same_cluster.append(line_num)
            line_num += 1    
        cluster.append(same_cluster)

    return cluster

def getKMeans(cluster,cluster_num,gene_id): 
    kmeans = [0,0,0,0,0,0,
                0,0,0,0,0,0,] #* 12-dimension
    object_num = len(cluster[cluster_num])

    for i in range(DIMENSION):
        
        tmp = 0
        for num in range(object_num):
            
            tmp += (float(gene_id[cluster[cluster_num][num]][i]))
            
        tmp /= float(object_num)
        
        kmeans[i] = float("{:.3f}".format(tmp))
        
    #*print("k-means : ", kmeans)
    return kmeans


def getDistance(object1, object2):
    dist = 0
    
    for dimension in range(DIMENSION):
        dist_dimesion = float(object1[dimension]) - float(object2[dimension])
        tmp = math.pow(dist_dimesion,2)
        dist += tmp
    dist = math.sqrt(dist)
    dist = "{:.3f}".format(dist)
    
    return float(dist)

def Reassign(obj, to_cluster_num,from_cluster_num,cluster):
    # * add to new cluster
    cluster[to_cluster_num].append(obj) # * number of object

    # * remove from old cluster
    #print("=====================")
    #print(cluster[from_cluster_num])
    #print(obj)
    cluster[from_cluster_num].remove(obj)
    #print(obj,"이", "cluster",to_cluster_num,"에 추가되었습니다.")
    #print(cluster[to_cluster_num])
    #print("------------")
    return cluster
    

def output_to_file(filename,cluster):
    file = open(filename, 'w')
    
    for i in range(k):
        file.write('{0}: '.format(len(cluster[i])))
        for v in cluster[i]:
            file.write(str(v)+" ")
        file.write("\n")
        

    file.close()
    print("Finished to print to output file : ", filename)

def SetCentroid(cluster,gene_id):
    cluster_num = 0
    centroid_list = [-1,-1,-1,-1,-1,
                    -1,-1,-1,-1,-1] #* the number of cluster is k:10, value:-1 for initializing.

    while(cluster_num < k):
        #print("DEBUG: Cluster ",cluster_num)
        kmeans = getKMeans(cluster,cluster_num=cluster_num,gene_id=gene_id)
        
        shortest_distance = getDistance(kmeans,gene_id[cluster[cluster_num][0]]) #!init value
        centroid = cluster[cluster_num][0] #!init value

        for i in range(len(cluster[cluster_num])):
            distance = getDistance(kmeans,gene_id[cluster[cluster_num][i]])
            
            if distance < shortest_distance:
                centroid = cluster[cluster_num][i]
                shortest_distance = distance

        centroid_list[cluster_num] = centroid
        #print("centroid : ",centroid) #? centroid means line number
        #print("shortest_distance :",shortest_distance)
        #print("===============================")
        cluster_num += 1
    print("centroid_list :",centroid_list)
    return centroid_list

def SetNewCluster(cluster,centroid_list,gene_id):
    isChanged = False
    for cluster_num in range(k):
        #print("Cluster",cluster_num,"=> ",cluster[cluster_num])
        idx = 0
        while idx < len(cluster[cluster_num]):
            
            
            if cluster[cluster_num][idx] == centroid_list[cluster_num]:
                #*print("This is the centroid.")
                idx += 1
                continue
            #*print("오브젝트 ",cluster[cluster_num][idx])
            
            distance_same = getDistance(gene_id[centroid_list[cluster_num]],gene_id[cluster[cluster_num][idx]])
            shortest_distance = distance_same #!init value
            for u in centroid_list:    
                distance_other = getDistance(gene_id[u],gene_id[cluster[cluster_num][idx]])

                if distance_other < shortest_distance:
                    to_cluster_number = centroid_list.index(u)
                    shortest_distance = distance_other
            #*print("shortest distance : ",shortest_distance)
            #*print("========================")
            
            if distance_same != shortest_distance:
                isChanged = True
                obj = cluster[cluster_num][idx]
                Reassign(obj,to_cluster_number,from_cluster_num = cluster_num,cluster=cluster)
                idx -= 1

            idx += 1

    return isChanged    

def Assigning(medoids,gene_id):
    isChagned = False
    len_data = len(gene_id)
    cluster = list()

    # * Make Clusters with Medoids
    for i in range(k):
        same_cluster = []
        same_cluster.append(medoids[i])
        cluster.append(same_cluster)
    #print(cluster)
    
    #* assigning each non-medoid to the nearest medoid.
    for i in range(len_data):
        
        if i in medoids:
            continue
        
        shortest_distance = getDistance(gene_id[medoids[0]],gene_id[i]) #! init value
        cluster_num = 0 #! init value
        for medoid in medoids:
            
            distance_other = getDistance(gene_id[medoid],gene_id[i])
            if distance_other < shortest_distance:
                shortest_distance = distance_other
                cluster_num = medoids.index(medoid)
                #print(i,"가 클러스터",cluster_num,"로 이동합니다.")
                isChagned = True
                
        
        cluster[cluster_num].append(i)

    return cluster,isChagned

def SetMedoids(medoids,cluster,gene_id):
    for _cluster_num in range(k):
        len_one_cluster = len(cluster[_cluster_num])
        sum_list = [0 for i in range(len_one_cluster)]
        for i in range(len_one_cluster):
            result = 0.0
            for j in range(len_one_cluster):
                if i==j:
                    continue
                result += getDistance(gene_id[cluster[_cluster_num][i]],gene_id[cluster[_cluster_num][j]]) 
            result = "{:.3f}".format(result)    
            sum_list[i] = float(result)
        #print(sum_list)
        #print(medoids)
        min_idx = sum_list.index(min(sum_list))
        #print("!",cluster[cluster_num][min_idx])
        medoids[_cluster_num] = cluster[_cluster_num][min_idx]
    return medoids

def Re_Assigning(cluster,medoids,gene_id):
    
    isChagned = False
    len_data = len(gene_id)

    
    #* assigning each non-medoid to the nearest medoid.
    for i in range(len_data):
        for clus_num in range(len(cluster)):
            
            if i in cluster[clus_num]:
            
                from_cluster_num = clus_num
        if i in medoids:
            continue
        
        shortest_distance = getDistance(gene_id[medoids[from_cluster_num]],gene_id[i]) #! init value
        to_cluster_num = from_cluster_num #! init value
        for medoid in medoids:
            
            distance_other = getDistance(gene_id[medoid],gene_id[i])
            if distance_other < shortest_distance:
                shortest_distance = distance_other
                to_cluster_num = medoids.index(medoid)
                if to_cluster_num != from_cluster_num:
                    isChagned = True
                
        if isChagned:
            if to_cluster_num == from_cluster_num:
                continue
            cluster = Reassign(i,to_cluster_num,from_cluster_num,cluster)
            
            
                
                
        
        

    return cluster,isChagned

def main():
    global DIMENSION,k
    DIMENSION = 12
    k = 10
    input_filename = 'assignment3_input.txt' #500
    #input_filename = 'test1.txt'
    output_filename = 'assignment3_output.txt'

    gene_id = getDataFromFile(input_filename)
    
    
    #print(gene_id)
    start_time = datetime.now()
    
    #? Step1-1. Calculate the distance between every pair of objects.
    len_data = len(gene_id)
    sum_list = [0 for i in range(len_data)]
    for i in range(len_data):
        result = 0.0
        for j in range(len_data):
            if i==j:
                continue
            result += getDistance(gene_id[i],gene_id[j]) #? Step1-2. Calculate the sum of distance for each object.

        result = "{:.3f}".format(result)    
        sum_list[i] = float(result)
    #print(sum_list)
    
    #? Step1-3. Select k objects having the smallest sum of distance as initial medoids.
    tmp_sum = sorted(sum_list)
    #print(tmp_sum[:k])
    #print(sum_list[:10])
    #print(sum_list.index(tmp_sum[0]))
    medoids = [-1 for i in range(k)]
    for i in range(k):
        medoids[i] = sum_list.index(tmp_sum[i])
    print("medoids : ",medoids)
    #? Step1-4. Obtain the initial clusters by assigning each non-medoid to the nearest medoid.
    cluster,isChanged = Assigning(medoids,gene_id)
    for i in range(k):
        print(cluster[i])
    
    isChanged = True
    count_for_debug = 0
    while isChanged : #?    3. Repeat steps 2-1 and 2-2 until the clusters do not change.
        #isChanged = False
        count_for_debug += 1
        print("[*]",count_for_debug,"번 반복하였습니다.")
        #? 1. For each cluster, calculate the sum of distance within the cluster for each object and select a new medoid having the smallest sum of distance.
        medoids = SetMedoids(medoids,cluster,gene_id)
        #*print(medoids)
        #?    2. Obtain the updated clusters by assigning each mon-medoid to the nearest medoid.
        cluster,isChanged = Re_Assigning(cluster,medoids,gene_id)
        
        
        
        
    
    print("\n")
    end_time = datetime.now()
    print("\n")
    size_ = 0
    for num in range(k):
        cluster[num].sort()
        print("cluster",num,"SIZE:",len(cluster[num]),cluster[num])
        size_ += len(cluster[num])
        print("\n")
    print("Total size : ", size_)
    output_to_file(output_filename,cluster)
    print("Time Elapsed : ", end_time - start_time,"microseconds")



if __name__ == '__main__':
    main()