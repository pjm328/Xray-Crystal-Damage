import os, pickle
import pandas as pd
import numpy as np

import img2cbf
reload(img2cbf)  # reload module for new changes(if any) to take effect
from img2cbf import img2cbf

def distl_wrapper(fake_cbf,distl_params):
	'''
	Takes a cbf_file and distl_params txt file
	Returns distl bash command to run
	e.g. distl_wrapper('fake_01_mos_0.10.cbf','distl_params.txt')
	'''
	with open(distl_params) as myfile:
		lines = myfile.readlines()
	lines = [x.strip() for x in lines]
	fixed_params = ['distl.signal_strength', fake_cbf, 'distl.verbose=True']
	return ' '.join(fixed_params+lines)


def ReadPeakIntensities(image_folder,header_contents,distl_params):
	'''
	Takes a folder of fake detector images with .img extension...
	and header_contents file describing simulation params...
	and distl_params file describing peak detection algorithm
	Returns a data frame with peak loc and intensity vs. frame no 
	1. Converts from img to cbf using img2cbf
	2. Finds peaks and calculates intensities using ReadPeakIntensities	
	3. Data frame is sorted by avg int of peaks over all frames 
	4. Writes the data frame to a pickle file for later use
	'''
	owd = os.getcwd()
	os.chdir(image_folder) # change wd for now, change back to owd at the end

	img_list = [x for x in os.listdir('.') if x.endswith('.img')]
	dic_peaks = {} # stores individual peak intensities
	eps = 5 # two peaks are considered same if within +/- eps in x and y

	for fnumber, fake_img in enumerate(img_list):
		# convert img to cbf 
		img2cbf(fake_img,header_contents,keep_original=True)
		fake_cbf = fake_img[0:-3]+'cbf'

		# run DISTL to find peaks using custom params in distl_params file
		# peak finding is highly sensitive to some particular params
		# more info: http://cci.lbl.gov/labelit/html/spotfinder.html
		bash_distl = distl_wrapper(fake_cbf,distl_params)		
		stdout = os.popen(bash_distl).read()
		
		# store individual peak intensities
		stdout_lines = stdout.splitlines()
		for line_num, line in enumerate(stdout_lines):
			if 'Total integrated signal=' in line:
				I = float(line.split()[7][7:])
				x = int(stdout_lines[line_num+2].split('x=')[1][0:4])
				y = int(stdout_lines[line_num+2].split('y=')[1][0:4])

				# check if (x,y) already in dic
				is_match = False
				for (x0,y0) in dic_peaks.keys():
					if x<=x0+eps and x>=x0-eps and y<=y0+eps and y>=y0-eps:
						is_match = True
						x_match = x0; y_match=y0
						break
				
				if is_match == False:
					dic_peaks[(x,y)] = [[fnumber],[I]]
				elif is_match == True:
					dic_peaks[(x_match,y_match)][0].append(fnumber)
					dic_peaks[(x_match,y_match)][1].append(I)


	# write peak intensities to a data frame for convenience
	df_peaks = pd.DataFrame( columns=['x','y'] + range(len(img_list)) )
	for peak in dic_peaks.keys():
		fnumbers,intensities = dic_peaks[peak]
		temp_dic = dict(zip(fnumbers,intensities))
		temp_dic['x'] = peak[0]
		temp_dic['y'] = peak[1]
		df_peaks = df_peaks.append(temp_dic,ignore_index=True)

	# calculate average intensity of each peak over all frames
	df_peaks['I-avg'] = np.nanmean(df_peaks.iloc[:,2::],axis=1)

	# sort peaks by the average intensity of peaks over all frames
	df_peaks.sort_values(by='I-avg', inplace=True, ascending=False)

	# save the data frame to a pickle file
	pickle.dump(df_peaks, open('df_peaks.pickle','wb'))

	os.chdir(owd)

	return df_peaks
