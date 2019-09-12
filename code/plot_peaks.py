import matplotlib.pyplot as plt
import pandas as pd
from get_peaks import *
from matplotlib.ticker import MultipleLocator

def plot_peaks(df_peaks,N,title,file_name):
    '''Plots N peaks from df_peaks with highest average I
    df_peaks - must be in form X Y I_1 I_2 ... I_avg
    N - number of peaks to plot
    title - title of plot
    file_name - name to save figure to (should be .svg)
    '''
    df_sorted = df_peaks.sort_values(by='I-avg',ascending=False)
    plt.close('all')
    fig, ax = plt.subplots()
    plt.ion()

    frames = df_sorted.columns[2:-1]
    for i in range(N):
        intensities = df_sorted.iloc[i][2:-1]
        plt.plot(list(frames),list(intensities))

    plt.xlabel('Frame number',fontsize='x-large')
    plt.ylabel('Peak intensity',fontsize='x-large')
    plt.title(title)
    plt.ylim(ymin = 0)
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    plt.show()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(file_name,format='svg')

    print("generated "+file_name)

data_folder = "/home/peter/Desktop/Summer_Research_2019/code/data/mseq6"
df_peaks = track_peaks(data_folder)
write_peaks("mseq6.csv",df_peaks)
plot_peaks(df_peaks,50,'Highest peaks with increasing mosaicity',
    'Highest_peaks_changing_mosaicity_6.svg')
