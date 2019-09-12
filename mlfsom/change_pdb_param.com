#! /bin/tcsh -f
set pdbfile = ""
set outfile = "modified.pdb"
set add_cell = ""
set add_bfactor = ""

foreach arg ($*)
	if (("$arg" =~ *.pdb)) then
		if (($pdbfile == "")) then
			if(-e "$arg") then
				set pdbfile = $arg
			else
	    	echo "WARNING: $arg does not exist."
			endif
		else
			set outfile = $arg
		endif
	endif
	if ("$arg" =~ add_cell=* && "$arg" =~ *[0-9]*) then
		set add_cell = `echo $arg | awk -F "=" '{print $NF+0}'`
	endif
	if ("$arg" =~ add_bfactor=* && "$arg" =~ *[0-9]*) then
		set add_bfactor = `echo $arg | awk -F "=" '{print $NF+0}'`
	endif
end

if (($pdbfile == "")||(($add_cell == "")&&($add_bfactor == ""))) then
	echo "ERROR: invalid arguments"
	exit 1
endif

cat $pdbfile | awk -v add_cell=$add_cell -v add_bfactor=$add_bfactor '/CRYST1/{A=$2+0;B=$3+0;C=$4+0;\
  newA=A+add_cell;newB=B+add_cell;newC=C+add_cell;\
	printf "CRYST1   %5.3f   %5.3f   %5.3f%s\n",newA,newB,newC,substr($0,34,71);next}\
	/^ATOM|^HETATM/{bfactor=$11+0;new_bfactor=bfactor+add_bfactor;\
	printf "%s%5.2f%s\n",substr($0,0,62),new_bfactor,substr($0,67,80);next} {print}' > $outfile
