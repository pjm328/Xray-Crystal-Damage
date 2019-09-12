#! /bin/bash
#

cp 1H87.pdb temp.pdb
#cp pdbfiles/4v1w.pd temp.pdb

# generate some fake structure factors from a Gd-soaked lysozyme structure
./ano_sfall.com energy=12660 refined.pdb 1.2A solvent_B=35

mv ideal_ano.mtz pristine.mtz

#grep -v " GD " refined.pdb >! damaged.pdb
#./ano_sfall.com energy=12660 damaged.pdb 1.2A solvent_B=35

#mv ideal_ano.mtz decayed.mtz

# maybe change some parameters
#vi mlfsom.params

start_mos=0.70
#bfactor_per_exposure=0.5
start_dose=.4
count=1
mos_inc=0.01
dose_inc=0.1
threads=1
#mos=0.2

mos=$start_mos
dose=$start_dose
i=0

rm /tmp/peter/*

#for ((i=1;i<=$count;i++)); do
for m in `cat failed_runs2.txt`; do
	mos=$m
	((i++))
	echo "mos=$mos runnning on $i"
	if [[ ! ($(($i%$threads)) -eq 0)]]; then
		#./change_pdb_param.com temp.pdb input$i.pdb cell_k=0.3 dose=$dose
		#./ano_sfall.com energy=12660 input$i.pdb 1.2A solvent_B=35
		#mv ideal_ano.mtz input$i.mtz
		#./mlfsom.com ~/Desktop/Summer_Research_2019/mlfsom/data/fake_cseq1-dose=${dose}-cellk=0.3_001.img input$i.pdb input$i.mtz mosaic=0.43 frames=1 osc=.01 &
		./mlfsom.com ~/Desktop/Summer_Research_2019/mlfsom/data/fake_mseq6-mos=${mos}_001.img pristine.mtz mosaic=$mos frames=1 osc=0.01 &

	else
		#./change_pdb_param.com temp.pdb input$i.pdb cell_k=0.3 dose=$dose
		#./ano_sfall.com energy=12660 input$i.pdb 1.2A solvent_B=35
		#mv ideal_ano.mtz input$i.mtz
		#./mlfsom.com ~/Desktop/Summer_Research_2019/mlfsom/data/fake_cseq1-dose=${dose}-cellk=0.3_001.img input$i.pdb input$i.mtz mosaic=0.43 frames=1 osc=.01 &
		./mlfsom.com ~/Desktop/Summer_Research_2019/mlfsom/data/fake_mseq6-mos=${mos}_001.img pristine.mtz mosaic=$mos frames=1 osc=0.01
		echo "joining threads..."
		wait
		echo "clearing temp files and moving results..."
		mv /tmp/peter/mlfsom*.XYI ~/Desktop/Summer_Research_2019/mlfsom/data/tempdata
		mv /tmp/peter/mlfsom*predin.txt ~/Desktop/Summer_Research_2019/mlfsom/data/tempdata
		for file in ~/Desktop/Summer_Research_2019/mlfsom/data/tempdata/mlfsom*predin.txt; do
		  tempnum=`echo $file | sed 's/[^0-9]*//g'`
		  mos=`cat $file | awk '/MOSAIC/{print $2}'`
		  rename "s/mlfsom_tempfile${tempnum}/mlfsom_tempfile_mseq6_mos=${mos}/" ~/Desktop/Summer_Research_2019/data/tempdata/*
		done
		for FILE in ~/Desktop/Summer_Research_2019/mlfsom/data/fake_mseq6-mos*.img; do
			drive upload -f $FILE -p 1b3_tNTErA5XcUsVxDke0tTRuYPezyZa_
		done
		for FILE in ~/Desktop/Summer_Research_2019/mlfsom/data/tempdata/*; do
			drive upload -f $FILE -p 1u7GvBtQWHEoAkxTAVPzwzjwb9aIQGD9k
		done
		rm ~/Desktop/Summer_Research_2019/mlfsom/data/fake_mseq6*.img
		#mv ~/Desktop/Summer_Research_2019/data/tempdata/* /media/peter/Lexar/Summer_Research_2019/data/tempdata
		rm ~/Desktop/Summer_Research_2019/mlfsom/fit2d_*
		rm /tmp/peter/*
		rm ~/Desktop/Summer_Research_2019/mlfsom/data/tempdata/*
	fi
	dose=`echo $dose $dose_inc | awk '{print $1+$2}'`
	mos=`echo $mos $mos_inc | awk '{print $1+$2}'`
done
