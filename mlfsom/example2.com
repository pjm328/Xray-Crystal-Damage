#! /bin/tcsh -f
#
#
#
#
cp 1H87.pdb refined.pdb

# generate some fake structure factors from a Gd-soaked lysozyme structure
./ano_sfall.com energy=12660 refined.pdb 1.2A solvent_B=35

mv ideal_ano.mtz pristine.mtz

#grep -v " GD " refined.pdb >! damaged.pdb
#./ano_sfall.com energy=12660 damaged.pdb 1.2A solvent_B=35

#mv ideal_ano.mtz decayed.mtz

# maybe change some parameters
#vi mlfsom.params

#./mlfsom.com ~/Desktop/Summer_Research_2019/mlfsom/data/fake_1_001.img frames=10
#rm /tmp/peter/*
./mlfsom.com ~/Desktop/Summer_Research_2019/mlfsom/data/fake_mos-4_001.img mosaic=.4 frames=1
#rm /tmp/peter/*
#./mlfsom.com ~/Desktop/Summer_Research_2019/mlfsom/data/fake_decay_001.img frames=10
rm /tmp/peter/*
#mv ~/Desktop/Summer_Research_2019/mlfsom/data/fake_decay_*.img /media/peter/Lexar/Summer_Research_2019/
rm ~/Desktop/Summer_Research_2019/mlfsom/fit2d_*

./mlfsom.com ~/Desktop/Summer_Research_2019/mlfsom/data/fake_mos-6_001.img mosaic=.6 frames=1
