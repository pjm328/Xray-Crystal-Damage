#! /bin/tcsh -f
set pdbfile = ""
set outfile = "modified.pdb"
set cell_k = ""
set dose=""

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
	if("$arg" =~ cell_k=* && "$arg" =~ *[0-9]*) then
		set cell_k = `echo $arg | awk -F "=" '{print $NF+0}'`
	endif
	if("$arg" =~ dose=* && "$arg" =~ *[0-9]*) then
		set dose = `echo $arg | awk -F "=" '{print $NF+0}'`
	endif
end

if (($pdbfile == "")||(($cell_k == "")&&($dose == ""))) then
	echo "ERROR: invalid arguments"
	exit 1
endif

cat $pdbfile | awk -v cell_k=$cell_k -v dose=$dose '/CRYST1/{A=$2+0;B=$3+0;C=$4+0;newA=A+cell_k*dose;newB=B+cell_k*dose;newC=C+cell_k*dose;\
	printf "CRYST1   %5.3f   %5.3f   %5.3f%s\n",newA,newB,newC,substr($0,34,71);next} {print}' > $outfile
