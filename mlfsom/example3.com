#! /bin/tcsh -f
#
#
#
#
cp 1H87.pdb refined.pdb

# generate some fake structure factors from a Gd-soaked lysozyme structure
./ano_sfall.com energy=12660 refined.pdb 1.2A solvent_B=35

mv ideal_ano.mtz pristine.mtz

# maybe change some parameters
#vi mlfsom.params

./mlfsom.com ./data/fake_2_001.img mosiac=3.2 frames=10

