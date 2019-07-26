#! /bin/tcsh -f
#
#
#
#alias rsh ssh -x
#
#  ssh-keygen
#  cp .ssh/id_rsa.pub .ssh/authorized_keys2
#
# shared directory that all computers can see
set sharedir = /tmp/${USER}/preds/
set localdir = /tmp/${USER}/
# 
mkdir -p ${sharedir}
mkdir -p ${localdir}

set filehost = `cd $sharedir ; df . | awk -F ":" '{print $1}' | tail -1`
if($#filehost != 1) set filehost = localhost
echo "files are being stored on $filehost"

if(! -e preds.txt) ln -sf ${sharedir}/preds.txt preds.txt


# generate the prediction file
set machines = `cat machines.txt`
if("$machines" == "") then
    set machines = "localhost"
endif


set Amatrix
set missets
set RESO = `echo "head" | mtzdump hklin pristine.mtz | awk '/Resolution Range/{getline;getline;print $6}'`
if("$RESO" == "") set RESO = 2.5
set CELL = `echo "head" | mtzdump hklin pristine.mtz | awk '/Cell Dimensions/{getline;getline;print}'`
if("$CELL" == "") then
    echo "error: cannot get cell from pristine.mtz"
    exit 9
endif
source mlfsom.params

foreach arg ( $* )
    set afternum = `echo $arg | awk -F "=" '$NF+0>0.0001{print $NF+0}'`
        if("$arg" =~ reso* && "$afternum" != "") then
	set temp = `echo $afternum | awk '$1>0.1 && $1<500{printf "%f", $1}'`
	if("$temp" != "") then
	    set RESO = "$temp"
	endif
    endif
end
echo "resolution $RESO"

set width = `echo $dim $pixel | awk '{print $1*$2}'`


if(! -e ref_predin.txt) then
    cat << EOF >! ref_predin.txt
# approximate a flat loop mu=4700 microns
ABSORB   70    0    0   4700
ABSORB  -70    0    0   4700
ABSORB    0  170    0   4700
ABSORB    0 -170    0   4700
ABSORB    0    0  200   4700
ABSORB    0    0 -600   4700
ABSORB    0  120  120   4700
ABSORB    0 -120  120   4700
ABSORB    0  150  -80   4700
ABSORB    0 -150  -80   4700
EOF
endif

awk 'NF>1' << EOF >! aux_predin.txt
RESOLUTION    $RESO
WAVELENGTH    $wavelength
DISPERSION    $dispersion
DIVERGENCE    $divergence
MOSAIC        $mosaic
SPEEDUP_PIX   $speedup_pix
SPEEDUP_PHI   $speedup_phi
MAX_SUB_SPOTS $max_subspots
POLARISATION  $polarization
DISTANCE      $distance
PIXEL         $pixel
DETSIZE       $width
BEAM          $beam_center
BEAMSIZE      $beam_size
XTALSIZE      $xtal_size
TWOTHETA      $twotheta
TILT          $detector_tilt
TWIST         $detector_twist
MISSETS       ${missets}
CELL          $CELL
MATRIX        $Amatrix
SPINDLE       $spindle_swing
AIR_MU        $mu_air
WINDOW_THICK  $window_thick
WINDOW_MU     $mu_window
DET_THICK     $det_thick
DET_MU        $mu_det
DET_MUEN      $muen_det
DET_MFP       $phosphor_mfp
EOF


rm -f busy

foreach machine ( $machines )
  rsh -n $machine "killall hkl2XYphi.awk ; rm ${localdir}/pred?" >& /dev/null
  sleep 0.1
end
echo "killed old jobs"
sleep 3
rm -f ${sharedir}/preds_cpu*.txt
echo "deleted old pred files"
echo "starting new jobs..."

set i = 0
set pwd = `pwd`
foreach machine ( $machines )

  @ i = ( $i + 1 )

  echo "$i $machine"
  echo "CPU $i $#machines" |\
  rsh $machine "cd $pwd ; cat ref_predin.txt aux_predin.txt - | ./hkl2XYphi.awk | dd of=${localdir}/pred$i " |& tee ${machine}_${i}.log &
  sleep 1

end
echo "waiting for jobs to finish..."
wait
echo "jobs are done."

# now do the copies one at a time to be nice to the file server
set i = 0
set pwd = `pwd`
foreach machine ( $machines )

@ i = ( $i + 1 )

echo "copying preds from $machine to $sharedir "
echo "$i $machine"
echo "" |\
 rsh $machine "mv ${localdir}/pred$i ${sharedir}/preds_cpu${i}_${machine}.txt" 

end




unset TOOBIG
unset BAD
unset TOOSMALL
set i = 0
foreach machine ( $machines )
    @ i = ( $i + 1 )
    set pred = ${sharedir}/preds_cpu${i}_${machine}.txt 
    echo -n "checking $pred "
    set test = `grep -l GOTHERE $pred | wc -l`
    if("$test") then
	echo "fubar"
	set BAD
	continue
    endif
    set test = `ls -l $pred | awk '{print ($5>=2.0^32)}'`
    if("$test") then
	echo "too big"
	set TOOBIG
	continue
    endif
    set test = `ls -l $pred | awk '{print ($5<10000)}'`
    if("$test") then
	echo "too small"
	set TOOSMALL
	continue
    endif
    echo "okay"
end
if($?TOOBIG || $?TOOSMALL || $?BAD) then
    echo "skipping"
    exit 9
endif



rm -f ${sharedir}/preds.txt
rsh -n $filehost "cat ${sharedir}/preds_cpu*.txt >  ${sharedir}/preds.txt"

echo "./preds.txt generated."

#rm -f ${sharedir}/preds_cpu*.txt

