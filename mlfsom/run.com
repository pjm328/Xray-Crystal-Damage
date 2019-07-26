#! /bin/bash
#

cp 1H87.pdb temp.pdb

# generate some fake structure factors from a Gd-soaked lysozyme structure
./ano_sfall.com energy=12660 refined.pdb 1.2A solvent_B=35

mv ideal_ano.mtz pristine.mtz

#grep -v " GD " refined.pdb >! damaged.pdb
#./ano_sfall.com energy=12660 damaged.pdb 1.2A solvent_B=35

#mv ideal_ano.mtz decayed.mtz

# maybe change some parameters
#vi mlfsom.params

start_mos=0.5
#bfactor_per_exposure=0.5
#start_exposure=3
count=1
mos_inc=0.05
#exposure_inc=2
threads=1
#mos=0.2

mos=$start_mos
exposure=$start_exposure

rm /tmp/peter/*

for ((i=1;i<=$count;i++)); do
	#addbfactor=`echo $exposure $bfactor_per_exposure | awk '{print $1*$2}'`
	echo "mos=$mos runnning on $i"
	if [[ ! ($(($i%$threads)) -eq 0)]]; then
		#./change_bfactor.com temp.pdb add=$addbfactor
		#mv modified_bfactor.pdb input$i.pdb
		#./ano_sfall.com energy=12660 input$i.pdb 1.2A solvent_B=35
		#mv ideal_ano.mtz input$i.mtz
		#./mlfsom.com ~/Desktop/Summer_Research_2019/mlfsom/data/fake_bseq1-${addbfactor}_001.img input$i.pdb input$i.mtz mosaic=$mos frames=1 osc=.01 &
		#inotifywait -re create /tmp/peter | while read path action file; do
		#	tempnum=$(echo $file | grep -o -E '[0-9]+')
		#	echo $tempnum $mos >> temp_file_key
		#done &
		./mlfsom.com ~/Desktop/Summer_Research_2019/mlfsom/data/fake_mseq4-mos=${mos}_001.img pristine.mtz mosaic=$mos frames=20 osc=1 &

	else
		#./change_bfactor.com temp.pdb add=$addbfactor
		#mv modified_bfactor.pdb input$i.pdb
		#./ano_sfall.com energy=12660 input$i.pdb 1.2A solvent_B=35
		#mv ideal_ano.mtz input$i.mtz
		#./mlfsom.com ~/Desktop/Summer_Research_2019/mlfsom/data/fake_bseq1-${addbfactor}_001.img input$i.pdb input$i.mtz mosaic=$mos frames=1 osc=.01
		#inotifywait -re create /tmp/peter | while read path action file; do
		#	tempnum=$(echo $file | grep -o -E '[0-9]+')
		#	echo $tempnum $mos >> temp_file_key
		#done &
		./mlfsom.com ~/Desktop/Summer_Research_2019/mlfsom/data/fake_mseq4-mos=${mos}_001.img pristine.mtz mosaic=$mos frames=20 osc=1
		echo "joining threads..."
		wait
		echo "clearing temp files and moving results..."
		mv /tmp/peter/mlfsom*.XYI ~/Desktop/Summer_Research_2019/data/tempdata
		mv /tmp/peter/mlfsom*predin.txt ~/Desktop/Summer_Research_2019/data/tempdata
		for file in ~/Desktop/Summer_Research_2019/data/tempdata/mlfsom*predin.txt; do
			tempnum=`echo $file | sed 's/[^0-9]*//g'`
			tempmos=`cat $file | awk '/MOSAIC/{print $2}'`
			rename "s/mlfsom_tempfile${tempnum}/mlfsom_tempfile_mos=${tempmos}/" ~/Desktop/Summer_Research_2019/data/tempdata/*
    		done
		rm ~/Desktop/Summer_Research_2019/mlsom/data/fake_mseq*.img
		#mv ~/Desktop/Summer_Research_2019/data/tempdata/* /media/peter/Lexar/Summer_Research_2019/data/tempdata
		rm ~/Desktop/Summer_Research_2019/mlfsom/fit2d_*
		rm /tmp/peter/*
	fi
	mos=`echo $mos $mos_inc | awk '{print $1+$2}'`
done

