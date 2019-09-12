import os, pickle
import pandas as pd
import numpy as np
from math import *

#Read in list of .XYI files and output a pd.DataFrame with
#format X Y I_1 I_2 ... I_avg

def get_files(XYI_folder):
    '''Takes a folder and returns list of files of .XYI format'''
    owd = os.getcwd()
    os.chdir(XYI_folder)
    path_name = os.getcwd()
    file_list = [path_name+'/'+x for x in os.listdir('.') if x.endswith('.XYI')]
    os.chdir(owd)
    return file_list

def get_peaks(peak_file_name,peak_dist=3):
    '''takes an .XYI file (and optionally peak clustering radius)
    returns a list of summed peaks where each peak is added together with its
    nearest neighbors within peak_dist
    returns peak_list - list of tuples in form (X,Y,Integrated_Intensity)
    '''
    #note: .XYI files have format - XDET YDET photons psi Rsize Tsize  h k l
    #parse peak file to get list of unclustered peaks
    peak_file = open(peak_file_name,'r')
    lines=peak_file.readlines()
    #peak_file.close()
    raw_peaks=[]
    for line in lines:
        split=line.rstrip().split(' ')
        raw_peaks.append((float(split[0]),float(split[1]),float(split[2])))

    #find clusters of peaks that are within peak_dist of eachother
    cluster_indicies=[0]
    for i in range(len(raw_peaks)-1):
        x1,y1=raw_peaks[i][0],raw_peaks[i][1]
        x2,y2=raw_peaks[i+1][0],raw_peaks[i+1][1]
        if sqrt((x1-x2)**2+(y1-y2)**2)>peak_dist:
            cluster_indicies.append(i+1)
    cluster_indicies.append(len(raw_peaks))

    #sum clusters together and return list of (X,Y,I_sum)
    peak_list=[]
    for i in range(len(cluster_indicies)-1):
        cluster=raw_peaks[cluster_indicies[i]:cluster_indicies[i+1]]
        weighed_x=0
        weighed_y=0
        total_I=0
        for peak in cluster:
            weighed_x+=peak[0]*peak[2]
            weighed_y+=peak[1]*peak[2]
            total_I+=peak[2]
        avg_x=weighed_x/total_I
        avg_y=weighed_y/total_I
        peak_list.append((avg_x,avg_y,total_I))
    print(len(peak_list))
    return peak_list

def track_peaks(XYI_folder,peak_dist=3,bin_size=20,image_height=2800,image_width=2800):
    '''Takes a folder of XYI files and outputs a pd.DataFrame with
    format X Y I_1 I_2 ... I_avg'''
    #get files
    file_names=get_files(XYI_folder)
    file_names.sort()
    #parse files and get lists of peaks
    frames=[]
    for file in file_names:
        frames.append(get_peaks(file,peak_dist=peak_dist))
        print(file+" read in")
    #split peaks into bins
    max_bin_i=image_width//bin_size+1
    max_bin_j=image_height//bin_size+1
    bins=[[[[] for f in range(len(frames))] for j in range(max_bin_j)] for i in range(max_bin_i)]
    #bin starting at x=i*bin_size,y=j*bin_size, for frame f, is given by
    #bins[i][j][f]
    print(image_height//bin_size+1,image_width//bin_size+1)
    for f in range(len(frames)):
        for p in range(len(frames[f])):
            peak=frames[f][p]
            x,y=peak[0],peak[1]
            i,j=int(x//bin_size),int(y//bin_size)
            if 0>i or i>max_bin_i or 0>j or j>max_bin_j:
                continue
            bins[i][j][f].append(p)
        print("frame "+str(f+1)+" put into bins")
    #track peaks through frames and put into DataFrame
    peak_used=[[False for peak in frame] for frame in frames]
    df_peaks=pd.DataFrame(columns=['x','y']+list(range(1,len(frames)+1)))
    for frame in range(len(frames)):
        for i in range(len(frames[frame])):
            if i%100==0:
                print(frame,i)
            if peak_used[frame][i]:
                continue
            peak_used[frame][i]=True
            peak=frames[frame][i]
            peak_dat=[peak[1],peak[0]]+[np.NaN]*frame+[peak[2]]
            for other_frame in range(frame+1,len(frames)):
                bin_i,bin_j=int(peak[0]//bin_size),int(peak[1]//bin_size)
                ind=[]
                for di in (-1,0,1):
                    for dj in (-1,0,1):
                        if bin_i+di<0 or bin_i+di>max_bin_i or bin_j+dj<0 or bin_j+dj>max_bin_j:
                            continue
                        ind+=bins[bin_i+di][bin_j+dj][other_frame]
                for j in ind:
                    if peak_used[other_frame][j]:
                        continue
                    other_peak=frames[other_frame][j]
                    if sqrt((other_peak[0]-peak[0])**2+(other_peak[1]-peak[1])**2)<peak_dist:
                        peak_dat.append(other_peak[2])
                        peak_used[other_frame][j]=True
                        break
            peak_dat+=[np.NaN]*(len(frames)+2-len(peak_dat))
            #if len(df_peaks)>54140: print(peak_dat)
            peak_dat.append(np.nanmean(peak_dat[2:len(peak_dat)]))
            temp_dic=dict(zip(['x','y']+list(range(1,len(frames)+1))+['I-avg'],peak_dat))
            df_peaks=df_peaks.append(temp_dic,ignore_index=True)
        print(str(frame+1)+" frame finished")
    return df_peaks

def write_peaks(file_name,df_peaks):
    '''Saves DataFrame containing peaks as csv in format:
    x,y,I_1,I_2,...,I_avg
    file_name - string specifying file name, should be a .csv
    df_peaks - DataFrame containing information about peaks, has columns in same
    format as csv
    '''
    file=open(file_name,'w')
    for i in range(len(df_peaks)):
        line=""
        params=df_peaks.iloc[i]
        for param in params:
            line+=str(param)+','
        if line!="":
            line=line[0:len(line)-1]+'\n'
        file.write(line)
    file.close()

def read_peaks(file_name):
    '''Reads a csv file and output DataFrame containing peak infomation
    file_name - string specifying file name, .csv file in format:
    x,y,I_1,I_2,...,I-avg
    '''
    file=open(file_name,'r')
    lines=file.readlines()
    file.close()
    frame_num=len(lines[0].rstrip().split(','))-3
    df_peaks=pd.DataFrame(columns=['x','y']+list(range(1,frame_num+1))+['I-avg'])
    i=0
    print(len(lines))
    for line in lines:
        i+=1
        if i%100==0:
            print(i)
        params=line.rstrip().split(',')
        peak_dat=[float(p) for p in params]
        temp_dic=dict(zip(['x','y']+list(range(1,frame_num+1))+['I-avg'],peak_dat))
        df_peaks=df_peaks.append(temp_dic,ignore_index=True)
    return df_peaks
