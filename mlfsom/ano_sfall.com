#! /bin/csh -f
#
#	Calculate data with anomalous differences using sfall, mapman and sftools		James Holton 5-14-14
#
#
#alias mapman /programs/o/rave/lx_mapman

# defaults
set pdbin   = ""
set sitepdb = ""
set fpp     = auto
set fp      = auto
set wave    = 1.000
set sigma   = 0.1

set outfile = ideal_ano.mtz

set SG    = ""
set reso  = ""
set solv_reso = 3
# rough "spacing" between grid points
set gridspacing = 0.38
set USER_GRID

set solvent_radius = 1.41
set solvent_scale  = 0.334
set solvent_B      = 50
set solsmooth_itrs = auto

set protein_scale  = 1.0
set protein_B      = 0

set partstruct_mtzs = ""

set site_B        = 0

set tempfile = /tmp/${USER}/tempfile$$

set hgen = 1

set logfile = ano_sfall_debug.log
rm -f $logfile
#set logfile = /dev/null

goto Setup
# read the command line
return_from_Setup:

# get cell dimensions
set CELL = `awk '/CRYST1/{print $2,$3,$4,$5,$6,$7}' $pdbin`
if("$SG" == "") then
    set pdbSG = `awk '/^CRYST/{print substr($0,56,12)}' $pdbin | head -1`
    if("$pdbSG" == "R 32") set pdbSG = "R 3 2"
    if("$pdbSG" == "P 21") set pdbSG = "P 1 21 1"
    if("$pdbSG" == "R 3 2" && $CELL[6] == 120.00) set pdbSG = "H 3 2"
    set SG = `awk -v pdbSG="$pdbSG" -F "[\047]" 'pdbSG==$2{print;exit}' ${CLIBD}/symop.lib | awk '{print $4}'`
    if("$SG" == R3 && $CELL[6] == 120.00) set SG = H3
endif
if("$SG" == "") set SG = P1


# get "site file" cell dimensions
set site_cell = `awk '$1 ~ /^CRYST/{printf "%.3f %.3f %.3f %.3f %.3f %.3f\n", $2,$3,$4,$5,$6,$7;exit}' $sitepdb`
if($#site_cell != 6 && $#CELL == 6) then
    echo "WARNING: no cell in $sitepdb"
    echo "         will use $CELL"
    set site_cell = ( $CELL )
endif
if($#site_cell == 6 && $#CELL != 6) then
    echo "WARNING: no cell in $pdbin"
    echo "         will use $site_cell"
    set CELL = ( $site_cell )
endif

# decide on a grid spacing
set default_res = `awk '/^REMARK   2 RESOLUTION/{print $4}' $pdbin`
if("$default_res" == "") then
    set minB = `awk '/^ATOM/ || /^HETATM/{print substr($0, 61, 6)+0}' $pdbin $sitepdb | sort -n | head -1`
    set default_res = `echo $minB | awk '$1>0{print 3*sqrt($1/80)}'`
endif
if("$default_res" == "") set default_res = 1.5

# reduce default resolution if it will overload sftools
set megapoints = `echo $CELL $default_res | awk '{printf "%d", 27*($1 * $2 * $3)/($NF*$NF*$NF)/900000}'`
if("$megapoints" > 100) then
    set default_res = `echo $CELL 100000000 | awk '{print ($1*$2*$3*27/$NF)^(1/3)}'`
endif
if("$reso" == "") then
    set reso = $default_res
endif

cat << EOF
pdbin $pdbin
sitepdb $sitepdb
reso $reso
cell $CELL
symm $SG
fp $fp
fpp $fpp
wave $wave
EOF

# cannonicalize the PDBs
echo "END" | cat $pdbin - >! ${tempfile}_model.pdb
echo "CELL $CELL" | pdbset xyzin ${tempfile}_model.pdb xyzout ${tempfile}.pdb >> $logfile
mv ${tempfile}.pdb ${tempfile}_model.pdb




# add hydrogens
if($hgen) then
    echo "adding hydrogens"
    cat ${tempfile}_model.pdb >! ${tempfile}.pdb
    echo "HYDROGENS SEPARATE" |\
    hgen xyzin ${tempfile}.pdb xyzout ${tempfile}_hydrogens.pdb  >> $logfile


    cat ${tempfile}.pdb |\
    awk '/^CRYST1/ || /^SCALE/{print}' |\
    cat >! ${tempfile}_model.pdb
    cat ${tempfile}.pdb ${tempfile}_hydrogens.pdb |\
    awk '/^ATOM/ || /^HETAT/' |\
    cat >> ${tempfile}_model.pdb
    echo "END" >> ${tempfile}_model.pdb
endif



# use SFALL to calculate absolute-scale electron density map
# pick a grid using SFALL
try_again:
echo "calculating electron density for the protein model"
sfall xyzin ${tempfile}_model.pdb mapout ${tempfile}_protein.map << EOF >> $logfile
MODE ATMMAP
CELL $CELL
SYMM $SG
resolution $reso
FORM NGAUSS 5
$USER_GRID
EOF
if($status) then
    # do something?
    set test = `grep "increase MEMSIZE" $logfile | wc -l`
    if("$test" == "1" && ! $?INCREASED_MEMSIZE) then
	echo "increasing MEMSIZE to maximum hugeness..."
	alias sfall sfall MEMSIZE 200000000
	set INCREASED_MEMSIZE
	goto try_again
    endif
    set BAD = "sfall failed to run"
    goto exit
endif

echo "go" | mapdump mapin ${tempfile}_protein.map | tee ${tempfile}maphead.txt | grep dens

set sfall_grid = `awk '/Grid sampling on x, y, z ../{print $8,$9,$10}' ${tempfile}maphead.txt  | tail -1`
echo "sfall chose grid: $sfall_grid"
set solv_grid = ( $sfall_grid )

# get the XYZ ordering
cat ${tempfile}maphead.txt |\
awk '/Fast, medium, slow / && / X /{print $(NF-2),$(NF-1),$NF}' |\
cat >! ${tempfile}sfallaxes.txt
set sfall_axes = `cat ${tempfile}sfallaxes.txt`
rm -f ${tempfile}sfallaxes.txt
echo "sfall_axes: $sfall_axes"

# get the grid on x,y,z
cat ${tempfile}maphead.txt |\
awk '/Start and stop points/{b["c"]=$(NF-5);e["c"]=$(NF-4);b["r"]=$(NF-3);e["r"]=$(NF-2);b["s"]=$(NF-1);e["s"]=$NF}\
     /Fast, medium, slow / && / X /{crs[$(NF-2)]="c";crs[$(NF-1)]="r";crs[$NF]="s";\
      print b[crs["X"]],e[crs["X"]],b[crs["Y"]],e[crs["Y"]],b[crs["Z"]],e[crs["Z"]];\
   }' |\
cat >! ${tempfile}xyzgrid.txt
set xyzgrid = `cat ${tempfile}xyzgrid.txt`
rm -f ${tempfile}xyzgrid.txt
echo "xyzgrid: $xyzgrid"

# get the grid on c,r,s
cat ${tempfile}maphead.txt |\
awk '/Start and stop points/{print $(NF-5),$(NF-4),$(NF-3),$(NF-2),$(NF-1),$NF}' |\
cat >! ${tempfile}crsgrid.txt
set crsgrid = `cat ${tempfile}crsgrid.txt`
rm -f ${tempfile}crsgrid.txt
echo "crsgrid: $crsgrid"



if($?USER_ANO) then
    if(! $?only_heavy) set only_heavy = 0

    # turn "site" PDB into anomalous atoms
    echo "END" | cat $sitepdb - |\
    awk -v switch=$only_heavy '! switch{print;next}\
        /^ATOM|^HETATM|^ANISO/ && ( substr($0, 13, 2) ~ / H| N| O| C/ || substr($0, 77, 2) ~ / H| N| O| C/){next} {print}' |\
    awk '/^ATOM/ || /^HETATM/{printf "ATOM%7d ANO  ANO A   1    %s\n",++n,substr($0,31,36);next}\
           /^ANISOU /{printf "ANISOU%5d ANO  ANO A   1    %s\n",n,substr($0,31,40);next} {print}' |\
    cat >! ${tempfile}_sites.pdb
    echo "CELL $CELL" | pdbset xyzin ${tempfile}_sites.pdb xyzout ${tempfile}.pdb >> $logfile
    mv ${tempfile}.pdb ${tempfile}_sites.pdb

    set elements = Ano
    echo "Ano $fp $fpp" >! ${tempfile}fp_fpp

else
    echo "using tabulated fp and fpp for atoms in $sitepdb"

    # make a simpler copy of the "site" pdb
    echo "END" | cat $sitepdb - |\
    awk '/^ATOM/ || /^HETATM/{print}' |\
    cat >! ${tempfile}_sites.pdb
    echo "CELL $CELL" | pdbset xyzin ${tempfile}_sites.pdb xyzout ${tempfile}.pdb >> $logfile
    mv ${tempfile}.pdb ${tempfile}_sites.pdb

    # see what elements are in here
    cat $sitepdb |\
    awk '/^ATOM/ || /^HETAT/{Ee=substr($0, 13, 2);\
	if(Ee=="  "){Ee=substr($0,13,4);gsub("[0-9]","",Ee)};\
	if(Ee~/^H/ && substr($0, 13, 3) !~ /^H[EFOG] /) Ee="H"; print Ee}' |\
    awk '{++sum[$1]} END{for(Ee in sum) print Ee}' >! ${tempfile}Ees
    set elements = `cat ${tempfile}Ees`
    rm -f ${tempfile}Ees

    # get the f' and f" values for each element at this wavelength
    rm -f ${tempfile}fp_fpp >& /dev/null
    echo "Ee  fp     fpp"
    foreach Ee ( $elements )
	set fp_fpp = `echo "$Ee $wave" | awk '{print "ATOM",$1;print "NWAVE 1",$2}' | crossec | awk 'NF==4 && $4+0>0{print $3,$4}'`

	echo "$Ee $fp_fpp" | tee -a ${tempfile}fp_fpp
    end

endif



if($?NO_SOLVENT) then
    goto skip_solvent
endif

# rudimentary bulk solvent

# make the all-or-nothing mask (with the electron density of H2O in e-/A^3)
# 1.00 g/cm^3 / 18.015 g/mol * 6.022e23 molecules/mol / (1e10 A/m / 100 cm/m)^3 * 10 e-/molecule = 0.334277
# note: sfall multiplies map by a basically random number!
#set solv_grid = `echo $CELL $solv_reso | awk '{print 2^int(log($1/$NF*3)/log(2)),2^int(log($2/$NF*3)/log(2)),2^int(log($3/$NF*3)/log(2))}'`
echo "making solvent mask: radius = $solvent_radius A  scale = $solvent_scale e-/A^3"
cavenv xyzin ${tempfile}_model.pdb mapout ${tempfile}raw_solvent.map << EOF >> $logfile
CELL $CELL
SYMM $SG
ENVSOLVENT
grid $sfall_grid
RADMAX $solvent_radius
EOF
set test = `ls -l ${tempfile}raw_solvent.map | awk '{print ($5>10000)}'`
if(! $test) then
    echo "WARNING: cavenv crashed! "
    echo "using sfall for solvent mask instead..."
    sfall xyzin ${tempfile}_model.pdb mapout ${tempfile}raw_solvent.map << EOF >> $logfile
    MODE ATMMAP SOLVMAP
    CELL $CELL
    SYMM $SG
    resolution $solv_reso
    grid $sfall_grid
    #H2OBG 0.33456
    vdwrad $solvent_radius
EOF
endif
# put map on the right scale
set bulk_level = `echo "" | mapdump mapin ${tempfile}raw_solvent.map | awk '/Maximum dens/{print $NF}'`
set scale = `echo "$solvent_scale $bulk_level" | awk '$2+0==0{print 0;exit} {print $1/$2}'`
echo "scale factor $scale 0" |\
mapmask mapin ${tempfile}raw_solvent.map mapout ${tempfile}.map >> $logfile
# and on same grid system as protein map
echo "axis $sfall_axes" |\
mapmask mapin ${tempfile}.map mapout ${tempfile}2.map >> $logfile
echo "xyzlim $xyzgrid" |\
mapmask mapin ${tempfile}2.map mapout sharp_solvent.map >> $logfile
rm -f ${tempfile}.map ${tempfile}2.map > /dev/null



echo "solvent map:"
echo "go" | mapdump mapin sharp_solvent.map | tee ${tempfile}maphead.txt | grep dens

# get the mean grid spacing
cat ${tempfile}maphead.txt |\
awk '/Grid sampling on x, y, z ./{as=$(NF-2);bs=$(NF-1);cs=$(NF-0)} \
     /Cell dimensions .../{a=$(NF-5);b=$(NF-4);c=$(NF-3);\
       print a/as,b/bs,c/cs,"mean =",(a/as + b/bs + c/cs) /3}' |\
cat >! ${tempfile}gridspacing.txt
set solvent_gridspacing = `cat ${tempfile}gridspacing.txt`
rm -f ${tempfile}gridspacing.txt
echo "solvent_gridspacing: $solvent_gridspacing"


# what B-factor will the smoothing effectively apply?
if("$solsmooth_itrs" == "auto") then
    # pick a reasonable number of iterations
    set solsmooth_itrs = `echo $solvent_B $solvent_gridspacing | awk '{print int($1/(2*($NF*3.1415*1.468)^2))}'`
    if("$solsmooth_itrs" == "") set solsmooth_itrs = 0
    if($solsmooth_itrs < 3) set solsmooth_itrs = 0
    echo "chose $solsmooth_itrs rounds of real-space solvent mask smoothing"
endif
set smoothB = `echo $solsmooth_itrs $solvent_gridspacing | awk '{print $1*2*($NF*3.1415*1.468)^2}'`
set rs_solvent_B = `echo $solvent_B $smoothB | awk '{print $1-$2}'`
echo "real-space smoothing will apply effective B of $smoothB"
echo "reciprocal-space B will be $rs_solvent_B"

# put solvent map on on XYZ for mapman
mapmask mapin1 sharp_solvent.map mapout ${tempfile}xyz_solvent.map << EOF >> $logfile
axis X Y Z
xyzlim ASU
EOF

# get the new grid so we can extend it a bit
echo "go" | mapdump mapin ${tempfile}xyz_solvent.map >! ${tempfile}maphead.txt
cat ${tempfile}maphead.txt |\
awk '/Start and stop points/{cs=$(NF-5);ce=$(NF-4);rs=$(NF-3);re=$(NF-2);ss=$(NF-1);se=$NF}\
     /Fast, medium, slow / && / X /{c=$(NF-2);r=$(NF-1);s=$NF;\
      print c,r,s,cs,ce,rs,re,ss,se;\
      print r,s,c,rs,re,ss,se,cs,ce;\
      print s,c,r,ss,se,cs,ce,rs,re;\
   }' |\
awk '/X Y Z/{print $4,$5,$6,$7,$8,$9}' |\
cat >! ${tempfile}xyzgrid.txt
set trimgrid = `cat ${tempfile}xyzgrid.txt`
rm -f ${tempfile}xyzgrid.txt
echo "trimgrid: $trimgrid"
set extgrid = `echo "$trimgrid" | awk '{print $(NF-5)-10,$(NF-4)+10,$(NF-3)-10,$(NF-2)+10,$(NF-1)-10,$NF+10}'`
echo "extgrid: $extgrid"

# expand the map a little so that smoothed edges will get trimmed off (mapman does not smooth with symmetry)
echo "expanding the solvent map slightly..."
echo "xyzlim $extgrid" |\
mapmask mapin1 ${tempfile}xyz_solvent.map mapout ${tempfile}ext.map >> $logfile

# determine how many times to apply the smoothing filter
echo "" >! ${tempfile}smooth.in
foreach itr ( `echo $solsmooth_itrs | awk '{for(i=1;i<=$1;++i) print i}'` )
echo "filter smooth map1 3" >> ${tempfile}smooth.in
end

# smooth it out
echo "smoothing solvent map in real space with mapman..."
rm -f ${tempfile}_smoothsolv.map >& /dev/null
setenv MAPSIZE `ls -l ${tempfile}ext.map | awk '{printf "%d", $5/3.5}'`
mapman -b mapsize $MAPSIZE << EOF >> $logfile
read map1 ${tempfile}ext.map CCP4
@${tempfile}smooth.in
write map1 ${tempfile}_smoothsolv.map CCP4
quit
y
EOF
if(! -e ${tempfile}_smoothsolv.map) then
    echo "WARNING: MAPMAN failed! "
    echo "using sharp-edged solvent mask"
    cp ${tempfile}ext.map ${tempfile}_smoothsolv.map
    set rs_solvent_B = $solvent_B
endif
echo "smoothed solvent map:"
echo "go" | mapdump mapin ${tempfile}_smoothsolv.map | grep dens
echo "go" | mapdump mapin ${tempfile}_smoothsolv.map | awk '/Start and stop points/{print "smoothsolv: ", $(NF-5),$(NF-4),$(NF-3),$(NF-2),$(NF-1),$NF}'


# trim off the garbled edges and put map axis ordering back into SFALL convention
echo "trimming off the extra edges..."
mapmask mapin1 ${tempfile}_smoothsolv.map mapout ${tempfile}trimmed.map << EOF >> $logfile
axis X Y Z
xyzlim $trimgrid
EOF
echo "go" | mapdump mapin ${tempfile}trimmed.map |\
 awk '/Start and stop points/{print "trimmed: ", $(NF-5),$(NF-4),$(NF-3),$(NF-2),$(NF-1),$NF}'

mapmask mapin1 ${tempfile}trimmed.map mapout ${tempfile}axes.map << EOF >> $logfile
axis $sfall_axes
EOF
#echo "go" | mapdump mapin ${tempfile}axes.map | awk '/Start and stop points/{print "axes: ", $(NF-5),$(NF-4),$(NF-3),$(NF-2),$(NF-1),$NF}'

mapmask MAPLIM ${tempfile}_protein.map mapin1 ${tempfile}axes.map mapout solvent_asu.map << EOF >> $logfile
# match appears to be broken
#xyzlim match
axis $sfall_axes
xyzlim $xyzgrid
EOF
rm -f ${tempfile}trimmed.map
rm -f ${tempfile}axes.map
echo "go" | mapdump mapin solvent_asu.map |\
  awk '/Start and stop points/{print "trim_crs: ", $(NF-5),$(NF-4),$(NF-3),$(NF-2),$(NF-1),$NF}'
echo "sfall_crs: $crsgrid"

# convert raw (sharp edge) solvent map to HKLs
echo "converting sharp solvent map into structure factors..."
sfall mapin sharp_solvent.map hklout ${tempfile}.mtz << EOF >> $logfile
MODE SFCALC MAPIN
CELL $CELL
SYMM $SG
resolution $reso
grid $sfall_grid
EOF
echo "stats nbin 20" | mtzdump hklin ${tempfile}.mtz |\
 awk '/PARTIAL FILE STAT/,""' |\
 awk '/FC/ && $1+0>0{print}\
  /Mean    Mean   Resolution/ && ! h1{h1=1;print} /Missing complete/ && ! h2{h2=1;print}' |\
 tee unsmooth_FC.txt

# convert smoothed solvent map to HKLs
echo "converting smooth solvent map into structure factors..."
sfall mapin solvent_asu.map hklout ${tempfile}.mtz << EOF >> $logfile
MODE SFCALC MAPIN
CELL $CELL
SYMM $SG
resolution $reso
grid $sfall_grid
EOF
echo "stats nbin 20" | mtzdump hklin ${tempfile}.mtz |\
awk '/PARTIAL FILE STAT/,""' |\
 awk '/FC/ && $1+0>0{print}\
  /Mean    Mean   Resolution/ && ! h1{h1=1;print} /Missing complete/ && ! h2{h2=1;print}' |\
 tee smooth_FC.txt

# apply a big B-factor
echo "applying B-factor of $rs_solvent_B (total B = $solvent_B) to the solvent structure factors"
cad hklin1 ${tempfile}.mtz hklout ${tempfile}_solvent.mtz << EOF >> $logfile
scale file 1 1 $rs_solvent_B
labin file 1 all
EOF


skip_solvent:



echo "calculating structure factors for the protein model"
sfall xyzin ${tempfile}_model.pdb hklout ${tempfile}.mtz << EOF >> $logfile
MODE SFCALC XYZIN
CELL $CELL
SYMM $SG
resolution $reso
FORM NGAUSS 5
grid $sfall_grid
EOF
echo "applying scale factor of $protein_scale and B=$protein_B to the protein structure factors"
cad hklin1 ${tempfile}.mtz hklout ${tempfile}_model.mtz << EOF >> $logfile
scale file 1 $protein_scale $protein_B
labin file 1 all
EOF

if($?NO_SOLVENT) then
    echo "setting all solvent density to zero..."
cad hklin1 ${tempfile}_model.mtz hklout ${tempfile}_solvent.mtz << EOF >> $logfile
scale file 1 0 0
labin file 1 all
EOF
    echo "making absolute_scale.map "
    cp ${tempfile}_protein.map absolute_scale.map

else
    echo "making absolute_scale.map "
    mapmask mapin1 ${tempfile}_protein.map mapin2 solvent_asu.map \
        mapout absolute_scale.map << EOF >> $logfile
    mode mapin1 mapin2
    maps add
EOF

endif




# make a "starting" partial structure file
cad hklin1 ${tempfile}_model.mtz hklout ${tempfile}_partstruct.mtz << EOF >> $logfile
scale file 1 0 0
labin file 1 all
EOF
# add any partial structures to this "blank slate"
foreach partstruct_mtz ( $partstruct_mtzs )
    echo "header" | mtzdump hklin $partstruct_mtz |\
    awk '/ Column Labels :/{getline;getline;for(i=1;i<=NF;++i)label[i]=$i}\
          / Column Types :/{getline;getline;for(i=1;i<=NF;++i){\
		if($i=="F" && F=="")F=label[i];if($i=="P" && P=="")P=label[i];}}\
	END{print F,P}' |\
    cat >! ${tempfile}F_P.txt
    set F   = `awk '{print $1}' ${tempfile}F_P.txt`
    set PHI = `awk '{print $2}' ${tempfile}F_P.txt`

    if("$F" == "") then
	set BAD = "unable to find an F in $partstruct_mtz"
	goto exit
    endif
    if("$PHI" == "") then
	# make one?
	echo "no phase in ${partstruct_mtz}, running refmac..."
	refmac5 hklin ${partstruct_mtz} xyzin $pdbin \
         hklout ${tempfile}refmacout.mtz xyzout ${tempfile}refmacout.pdb << EOF | tee -a $logfile | grep "R factor"
	refi type rigid
	rigid ncyc 1
EOF
        echo "adding refmac-scaled $F from $partstruct_mtz with PHIC_ALL to partial structure"
        cad hklin1 ${tempfile}refmacout.mtz hklout ${tempfile}_part.mtz << EOF >> $logfile
        labin file 1 E1=FP E2=PHIC_ALL
EOF
	if($status) then
	    set BAD = "unable to obtain a phase for $partstruct_mtz"
	    goto exit
	endif
    else
        echo "adding $F $PHI in $partstruct_mtz to partial structure"
        cad hklin1 $partstruct_mtz hklout ${tempfile}_part.mtz << EOF >> $logfile
        labin file 1 E1=$F E2=$PHI
EOF
    endif

    sftools << EOF >> $logfile
read ${tempfile}_partstruct.mtz
read ${tempfile}_part.mtz
set labels
Fpart
PHIpart
FC
PHIC
calc ( COL Fpart PHIpart ) = ( COL Fpart PHIpart ) ( COL FC PHIC ) +
write ${tempfile}_new.mtz col Fpart PHIpart
y
stop
EOF
    mv ${tempfile}_new.mtz ${tempfile}_partstruct.mtz
end


# use SFALL to calculate "site" Fs
rm -f ${tempfile}_sites.mtz
foreach Ee ( $elements )

    # get the names right
    set pdbEe = $Ee
    if("$pdbEe" == "Ano") set pdbEe = AN

    # extract fp and fpp from the records
    set fp = `awk -v Ee="$Ee" '$1==Ee{print $2}' ${tempfile}fp_fpp`
    set fpp = `awk -v Ee="$Ee" '$1==Ee{print $3}' ${tempfile}fp_fpp`
    set test = `echo $fp $fpp | awk '{print ($1*$1>0 || $2*$2>0)}'`
    if(! $test) then
	echo "WARNING: unknown element: $Ee"
	echo "fp and fpp unknown.  Setting to zero"
	continue
    endif

    # extract appropriate atoms from the "site" file
    cat ${tempfile}_sites.pdb |\
    awk -v Ee="$pdbEe" '{ele=substr($0,13,2);\
	if(ele=="  "){ele=substr($0,13,4)};gsub("[0-9 ]","",ele);\
	if(ele~/^H/ && substr($0,13,3) !~ /^H[EFOG] /) ele=="H"}\
	 /^ATOM|^HETAT|^ANISOU/ && Ee!=ele{next} {print}' |\
    awk '/^ATOM/ || /^HETATM/{printf "ATOM%7d ANO  ANO A   1    %s\n",++n,substr($0,31,36);next}\
           /^ANISOU /{printf "ANISOU%5d ANO  ANO A   1    %s\n",n,substr($0,31,40);next} {print}' |\
    cat >! ${tempfile}_element.pdb
    echo "CELL $CELL" | pdbset xyzin ${tempfile}_element.pdb xyzout ${tempfile}.pdb >> $logfile
    mv ${tempfile}.pdb ${tempfile}_element.pdb

    set atoms = `awk '/^ATOM/' ${tempfile}_element.pdb | wc -l`
    echo "calculating anomalous structure factor contribution from $atoms $Ee : fpp=$fpp fp=$fp ..."

    sfall xyzin ${tempfile}_element.pdb hklout ${tempfile}_element.mtz << EOF >> $logfile
MODE SFCALC XYZIN
CELL $CELL
SYMM $SG
resolution $reso
GRID $sfall_grid
FORMFACTOR NGAUSS 5 Ano
EOF

    # use cad to sanitize file and perhaps apply a big B-factor to the sites
    mv ${tempfile}_element.mtz ${tempfile}.mtz >& /dev/null
    cad hklin1 ${tempfile}.mtz hklout ${tempfile}_element.mtz << EOF >> $logfile
scale file 1 1.0 $site_B
labin file 1 all
EOF

    if(! -e ${tempfile}_sites.mtz) then
	echo "starting new anomalous file"
	sftools << EOF >> $logfile
read ${tempfile}_element.mtz
set labels
FC
PHIC
calc ( COL Ffp PHIfp ) = ( COL FC PHIC ) [ $fp 0 ] *
calc ( COL Ffpp PHIfpp ) = ( COL FC PHIC ) [ 0 $fpp ] *
write ${tempfile}_sites.mtz col Ffp PHIfp Ffpp PHIfpp
y
stop
EOF
    else
	echo "adding to total..."
	sftools << EOF >> $logfile
read ${tempfile}_sites.mtz
read ${tempfile}_element.mtz
set labels
Fsum_fp
PHIsum_fp
Fsum_fpp
PHIsum_fpp
FC
PHIC
calc ( COL Ffp PHIfp    ) = ( COL FC PHIC ) [ $fp 0 ] *
calc ( COL Ffpp PHIfpp ) = ( COL FC PHIC ) [ 0 $fpp ] *
calc ( COL Ffp PHIfp    ) = ( COL Ffp PHIfp ) ( COL Fsum_fp PHIsum_fp ) +
calc ( COL Ffpp PHIfpp ) = ( COL Ffpp PHIfpp ) ( COL Fsum_fpp PHIsum_fpp ) +
write ${tempfile}_new.mtz col Ffp PHIfp Ffpp PHIfpp
y
stop
EOF
        mv ${tempfile}_new.mtz ${tempfile}_sites.mtz
    endif

end


# add these together properly
sftools << EOF >> $logfile
read ${tempfile}_model.mtz
read ${tempfile}_sites.mtz
read ${tempfile}_solvent.mtz
read ${tempfile}_partstruct.mtz
set labels
FP
PHIP
FHp
PHIHp
FHpp
PHIHpp
Fsolv
PHIsolv
Fpart
PHIpart
# add the solvent contribution to the protien
calc ( COL FP  PHIP ) = ( COL FP PHIP ) ( COL Fsolv PHIsolv ) +
# add the partial structure contribution to the protien
calc ( COL FP  PHIP ) = ( COL FP PHIP ) ( COL Fpart PHIpart ) +
# add dispersive component to protein to get FPH
calc ( COL FPH PHIPH ) = ( COL FP PHIP ) ( COL FHp PHIHp ) +
# add/subtract anomalous component to each Friedel mate
calc ( COL Fplus  PHIplus  ) = ( COL FPH PHIPH ) ( COL FHpp PHIHpp ) +
calc ( COL Fminus PHIminus ) = ( COL FPH PHIPH ) ( COL FHpp PHIHpp ) -
# compute difference to get dano
calc D COL DANO = COL Fplus COL Fminus -
# compute average to get Fmean
calc F COL Fsum = COL Fplus COL Fminus +
calc F COL Fmean = COL Fsum 2 /
delete Fsum
# compute total FH contribution?
calc ( COL FH PHIH ) = ( COL FHp PHIHp ) ( COL FHpp PHIHpp ) +
# assign phony errors
calc Q COL SIGFP = $sigma
calc Q COL SIGFPH = $sigma
calc Q COL SIGDANO = COL SIGFPH 2 *
write ${tempfile}_combo.mtz col FP SIGFP PHIP FH PHIH FPH SIGFPH DANO SIGDANO Fplus Fminus Fmean FHp PHIHp FHpp PHIHpp
y
stop
EOF

# sanitize the output with CAD
cad hklin1 ${tempfile}_combo.mtz hklout $outfile << EOF >> $logfile
LABIN file 1 E1=FP E2=SIGFP E3=PHIP E4=FH E5=PHIH E6=FPH E7=SIGFPH \
             E8=DANO E9=SIGDANO E10=Fplus E11=Fminus E12=Fmean \
             E13=FHp E14=PHIHp E15=FHpp E16=PHIHpp
CTYP  file 1 E1=F  E2=Q     E3=P    E4=F  E5=P    E6=F   E7=Q      \
             E8=D    E9=Q       E10=G     E11=G      E12=F  \
             E13=F E14=P E15=F E16=P
EOF

# check to see if this worked
set test = `echo "head" | mtzdump hklin $outfile | awk '/Number of Columns/{print $NF}'`
if("$test" != "19") then
    set BAD = "something went wrong.  Output mtz is corrupt."
endif

exit:
if($?BAD) then
    echo "ERROR: $BAD"
    exit 9
endif

# calculate Bijovet ratio
echo "" |\
mtzdump hklin $outfile |\
awk 'NF>3 && $NF=="DANO" && $(NF-1)=="D"{DANO=substr($0,48)+0} \
     NF>3 && $NF=="Fmean" && $(NF-1)i=="F"{Fmean=substr($0,48)+0}\
     END{print "Bijvoet ratio:",DANO"/"Fmean,"=",100*DANO/Fmean"%"}'

cat << EOF
rendered model in $pdbin
at ${reso}A with cell $CELL in $SG
with anomalous positions from $sitepdb
sigmas abritrarily assiged to $sigma
EOF
if($?USER_ANO) then
echo "fpp = $fpp"
echo "fp  = $fp"
else
echo "fp and fpp at wave = $wave"
endif

cleanup:
# clean up
#exit
rm -f ${tempfile}.map ${tempfile}.mtz >& /dev/null
rm -f ${tempfile}_model.pdb ${tempfile}_sites.pdb >& /dev/null
rm -f ${tempfile}_model.mtz ${tempfile}_sites.mtz >& /dev/null
rm -f ${tempfile}_solvent.mtz ${tempfile}_smoothsolv.mtz >& /dev/null
rm -f ${tempfile}_wrong.mtz >& /dev/null
rm -f ${tempfile}_sites.mlphare ${tempfile}_mlphareFH.mtz >& /dev/null
rm -f ${tempfile}_added.mtz ${tempfile}_addup.mtz >& /dev/null
rm -f ${tempfile}.pdb >& /dev/null
rm -f ${tempfile}_model.pdb >& /dev/null
rm -f ${tempfile}_protein.map >& /dev/null
rm -f ${tempfile}smooth.in >& /dev/null
rm -f ${tempfile}_smoothsolv.map >& /dev/null
rm -f ${tempfile}ext.map >& /dev/null
rm -f ${tempfile}fp_fpp >& /dev/null
rm -f ${tempfile}xyz_solvent.map >& /dev/null
rm -f ${tempfile}_hydrogens.pdb >& /dev/null
rm -f ${tempfile}raw_solvent.map >& /dev/null
rm -f ${tempfile}_element.pdb >& /dev/null
#rm -f solvent_asu.map >& /dev/null
#rm -f sharp_solvent.map >& /dev/null
rm -f ${tempfile}maphead.txt >& /dev/null
rm -f ${tempfile}_combo.mtz >& /dev/null
rm -f ${tempfile}_partstruct.mtz >& /dev/null
rm -f ${tempfile}_element.mtz >& /dev/null
rm -f ${tempfile}_sites.pdb >& /dev/null

exit


Setup:

set i = 0
foreach arg ( $* )
    @ i = ( $i + 1 )
    if( "$arg" =~ *.pdb || "$arg" =~ *.brk ) then
	if(-e "$arg") then
	    if("$pdbin" == "") then
		set pdbin = $arg
	    else
		set sitepdb = $arg
	    endif
	else
	    echo "ERROR: $arg does not exist"
	endif
    endif
    if( "$arg" =~ *.mtz ) then
	if(-e "$arg") then
	    set partstruct_mtzs = ( $partstruct_mtzs $arg )
	else
	    echo "ERROR: $arg does not exist"
	endif
    endif
    if( "$arg" =~ *A ) then
	set test = `echo $arg | awk '$1+0>0.1{print $1+0}'`
    	if("$test" != "") set reso = $test
    endif
    if( "$arg" =~ reso*=* ) then
	set test = `echo $arg | awk -F "[=]" '$2>0.1{print $2+0}'`
    	if("$test" != "") set reso = $test
    endif
    if( "$arg" == "-grid" ) then
	@ j = ( $i + 1 )
	@ k = ( $i + 2 )
	@ l = ( $i + 3 )
	if( $l > $#argv ) then
	    echo "WARNING: bad grid on command line"
	    continue
	endif
	set test = `echo "$argv[$j] $argv[$k] $argv[$l]" | awk '$1*$2*$3>0{print}'`
    	if("$test" != "") then
	    echo "user-specified grid $test"
	    set USER_GRID = "GRID $test"
	else
	    echo "WARNING: ignoring grid "$argv[$j] $argv[$k] $argv[$l]"
	endif
    endif
    if( "$arg" =~ wave=* ) then
	set USER_WAVE
	set test = `echo $arg | awk -F "[=]" '$2>0.05{print $2+0}'`
    	if("$test" != "") set wave = $test
    endif
    if( "$arg" =~ energy=* ) then
        set USER_WAVE
        set test = `echo $arg | awk -F "[=]" '$2<200 && $2>1{print 12.3984245/$2}'`
        if("$test" != "") set wave = $test
        set test = `echo $arg | awk -F "[=]" '$2>1000{print 12398.4245/$2}'`
        if("$test" != "") set wave = $test
    endif
    if( "$arg" =~ fp=* ) then
	set USER_ANO
	set test = `echo $arg | awk -F "[=]" '$2*$2>0{print $2+0}'`
    	if("$test" != "") set fp = $test
    endif
    if( "$arg" =~ fpp=* ) then
	set USER_ANO
	set test = `echo $arg | awk -F "[=]" '$2+0>0{print $2+0}'`
    	if("$test" != "") set fpp = $test
    endif
    if( "$arg" =~ sigma=* ) then
	set test = `echo $arg | awk -F "[=]" '$2+0>0{print $2+0}'`
    	if("$test" != "") set sigma = $test
    endif
    if( "$arg" =~ site_B=* ) then
	set test = `echo $arg | awk -F "[=]" '$2+0!=0{print $2+0}'`
    	if("$test" != "") set site_B = $test
    endif
    if( "$arg" =~ solvent_B=* || "$arg" =~ b_sol=* || "$arg" =~ bsol=*  || "$arg" =~ Bsol=* ) then
	set test = `echo $arg | awk -F "[=]" '$2~/^[0-9]/{print $2+0}'`
    	if("$test" != "") set solvent_B = $test
    endif
    if( "$arg" =~ solv_reso=* ) then
	set test = `echo $arg | awk -F "[=]" '$2~/^[0-9]/{print $2+0}'`
    	if("$test" != "") set solv_reso = $test
    endif
    if( "$arg" =~ solsmooth_itrs=* ) then
	set test = `echo $arg | awk -F "[=]" '$2+0>=0{print $2+0}'`
    	if("$test" != "") set solsmooth_itrs = $test
    endif
    if( "$arg" =~ solvent_scale=* ) then
	set test = `echo $arg | awk -F "[=]" '$2+0>=0{print $2+0}'`
    	if("$test" != "") set solvent_scale = $test
    endif
    if( "$arg" =~ k_sol=* || "$arg" =~ ksol=*  ) then
	set test = `echo $arg | awk -F "[=]" '$2+0>=0{print $2+0}'`
    	if("$test" != "") then
	    set solvent_scale = $test
	endif
    endif
    if( "$arg" =~ solvent_radius=* ) then
	set test = `echo $arg | awk -F "[=]" '$2+0>0{print $2+0}'`
    	if("$test" != "") set solvent_radius = $test
    endif
    if( "$arg" =~ scale=* ) then
	set test = `echo $arg | awk -F "[=]" '$2+0>=0{print $2+0}'`
    	if("$test" != "") set protein_scale = $test
    endif
    if( "$arg" =~ B=* ) then
	set test = `echo $arg | awk -F "[=]" '$2+0!=0{print $2+0}'`
    	if("$test" != "") set protein_B = $test
    endif
    if( "$arg" == nohydrogen ) then
	set hgen = 0
    endif
    if( "$arg" == nosolvent ) then
	set NO_SOLVENT
    endif
    if( "$arg" =~ [PpCcIiFfRrHh][1-6]* ) then
	set SG = "$arg"
    endif
end

if(! -e "$pdbin") then
    cat << EOF
usage: $0 model.pdb sites.pdb  P212121 2A fpp=4 fp=-8 sigma=1

where:
model.pdb is a brookhaven file of the protein structure
sites.pdb is a brookhaven file containing position/occupancy/B of anomalous scatterers

optional:
P212121   is the desired space group
2A        render to a desired resolution
wave=0.98 calculate the anomalous f' and f'' contribution at 0.98 A
fpp=4     sets the anomalous f'' contribution to 4 electrons
fp=-8     sets the dispersive f' contribution to -8 electrons
grid=0.3  use a particular grid spacing in sfall
sigma=1   gives a value of 1 to the sigma of each column


EOF
exit 9
endif

# default to same PDB for "sites"
if("$sitepdb" == "") then
    set sitepdb = $pdbin
    # ignore H,C,N,O atoms
    set only_heavy = 1
endif



goto return_from_Setup
