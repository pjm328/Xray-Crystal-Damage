#! /bin/tcsh -f
#
#
#
#

# generate some fake data from a Gd-soaked lysozyme structure
./ano_sfall.com energy=12660 1H87.pdb

# take ARP/wARP refined model
cat ~jamesh/projects/example_data_sets/ALS831/lyso1/refmac/files/molrep.pdb |\
grep -v " 0.00 " |\
awk '{gsub("DUM  DUM"," O   HOH"); print}' |\
cat > ! refined.pdb

# generate some fake data from a Gd-soaked lysozyme structure
./ano_sfall.com energy=12660 refined.pdb 1.2A solvent_B=35

mv ideal_ano.mtz pristine.mtz

grep -v " GD " refined.pdb >! damaged.pdb
./ano_sfall.com energy=12660 damaged.pdb 1.2A solvent_B=35

mv ideal_ano.mtz decayed.mtz

vi mlfsom.params

./mlfsom.com frames=100

