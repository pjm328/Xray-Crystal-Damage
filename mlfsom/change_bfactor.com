#! /bin/tcsh -f
set pdbfile = ""
set bfactormult = ""
foreach arg ($*)
	if ("$arg" =~ *.pdb) then
		if(-e "$arg") then
			set pdbfile = $arg
		else
	    		echo "WARNING: $arg does not exist."
		endif
	endif
	if("$arg" =~ mult=* && "$arg" =~ *[0-9]*) then
		set bfactormult = `echo $arg | awk -F "=" '{print $NF+0}'`
	endif
	if("$arg" =~ add=* && "$arg" =~ *[0-9]*) then
		set bfactoradd = `echo $arg | awk -F "=" '{print $NF+0}'`
	endif
end

if (($pdbfile == "")||(($bfactormult == "")&&($bfactoradd == ""))) then
	echo "ERROR: invalid arguments"
	exit 1
endif

if ($bfactormult == "") then
	set bfactormult = "1"
endif
if ($bfactoradd == "") then
	set bfactoradd = "0"
endif

cat $pdbfile | awk -v mult=$bfactormult -v add=$bfactoradd '/^ATOM/ || /^HETAT/{B=substr($0, 61, 6)+0;newB=B*mult+add;\
	printf "%s%6.2f%s\n",substr($0, 0, 61),newB,substr($0,67,12);next} {print}' > modified_bfactor.pdb


