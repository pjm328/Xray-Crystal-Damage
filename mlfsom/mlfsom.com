#! /bin/tcsh -f
#
#	convert a merged .mtz file into ADSC formatted images	- James Holton 1-31-16
#
#
if(-e mlfsom.sourceme) then
    echo "reading mlfsom.sourceme"
    source mlfsom.sourceme
endif

set outprefix    = /data/peter/fakedata/lyso_from_pdb/fake_1

set tempfile = ${CCP4_SCR}/mlfsom_tempfile$$
#set tempfile = ${CCP4_SCR}/mlfsom_tempfile
#set tempfile = /dev/shm/mlfsom_tempfile
#set tempfile = mlfsom_tempfile

unset CACHE

set path = ( /programs/FIT2D $path `dirname $0` )

set gnuplotversion = `gnuplot --version | awk '{print int($2)}'`
if($gnuplotversion == "") then
    echo "ERROR: need gnuplot installed to render backgrounds."
    exit 9
endif
set gnuploterminal = 'set terminal table'
if($gnuplotversion >= 5) set gnuploterminal = ""

set mtzfile         = pristine.mtz
set decayedmtz      = decayed.mtz
set decay_halfdose  = 5MGy
set decay_doseratio = 2000ph/um^2/Gy
set decay_AperMGray = 0.1386
set starting_dose   = 0MGy

set dimx           = 3072pix
set dimy           = 3072pix
set bitsperpixel   = 16
set pixel          = 0.1024mm
set EOgain         = 7.3electrons/photon
set ampgain        = 4electron/adu
# expected values for ADSC: Gd2O2S:Tb phospor and Al/plastic window
set window_thick   = 135.4
set mu_window      = 2151
set det_thick      = 32
set mu_det         = 18.6
set muen_det       = 20.3

# mean-free path of a visible photon in the phosphor
set phosphor_mfp   = 100

set extra_A        = 0
set extra_B        = 0
set peak_dqe       = 1.0
set vignette       = 0.4
set readnoise      = 11.5electrons/pixel
set ripplenoise    = 1%
set ripplenoise_size = 10
set adcoffset      = 40adu

set psf            = 0.06mm
set spotwidth      = 1.0pix
set spot_type      = fiber
set fiber_spread   = 30um
set fiber_length_ratio = 1.0865
# minimum width of plotted Gaussians in pixels
set speedup_pix    = 1.5
# minimum width of plotted Gaussians in mosaic widths
set speedup_phi    = 1


set beam_size    = 100um
set wavelength   = 1.115879A
set dispersion   = 0.00025
set divergence   = ( 0.10 0.02 )
set polarization = 0.9

set flux         = 2e+11ph/s
set exposure     = 1s

# rms variation of intensity for 1s exposure
set beam_flicker   = 0.1%/sqrtHz
set shutter_jitter = 0.5774ms
# dead time due to read-out in so-called-shutterless mode. positive is a gap, negative makes overlap (non physical?)
set shutter_gap    = 0 ms

set distance    = 100mm
set twotheta    = 0deg
set beam_center = ( 108.66 102.246 )
set beam_center = ( 157 157 )
set detector_tilt  = 0.3deg
set detector_twist = 0.1deg
set beamstop_size  = 0.3mm
set beamstop_dist  = 35mm

set xtal_size    = 100um
set missets     = ( 12.345 67.89 101.112 )
set missets     = ( 0 0 0 )
set mosaic      = 0.5deg
set rocking_curve = tanh
set phi0        = 0deg
set osc         = 1deg
set spindle_swing = 0deg
set range       = 1


# background scatter contributions: external files must be available!
set nonxtal_materials   = ( air     water Paratone-N nanoice iceIh )
# thickness crossed by the beam
set nonxtal_thickness   = ( 37e3    150   150          0       0   )
# density in g/cm^3
set nonxtal_densities   = ( 1.88e-3 1.0   0.8        1.0     0.92  )
# formula weight in amu (corresponding to "electron" units in file)
set nonxtal_molweight   = ( 28      18    12         18      18    )


# crystal incoherent elastic scattering?
# this depends on the B factor of the atoms in the PDB
set pdbfile = refined.pdb
set decay_A = 0

# crystal Compton scattering?
# this depends on the atomic composition of the PDB and the wavelength

# crystal fluorescence?
# this depends on the atomic composition of the PDB, the wavelength, and the fluorescent yield


# arbitrary... need to work on physical origin
set SAXS_intensity       = 6000adu/s
set icering_intensity = 0adu
set icering_width     = 1.5

# dead zones of a tiled CCD detector
set windowpanes_x   = 2
set windowpanes_y   = 2
set windowpane_width_x = 3pix
set windowpane_width_y = 3pix
set windowpane_value = 0

# debugging
#set frames = 3
#set spot_type    = pinpoints
#set spot_type    = gauss
#set spot_type    = lg
#set mosaic      = 0.01
#set spotwidth    = 1
#set NO_NOISE
#set NO_BACKGROUND
#set NO_SPOTS
#set icering_intensity = 0
#set SAXS_intensity = 0

set start_time = `msdate.com | awk '{print $NF}'`


# calculate (electrons/photon) / (electrons/ADU) =  ADU/photon
set gain = `echo $EOgain $ampgain | awk '{print $1/$2}'`


goto Setup
Return_from_Setup:

if ($#argv == 0) then
    cat << EOF

usage: mlfsom.com file_1_001.img frames=100 [1s] [200mm] [12660eV] [param=value]

where:
file_0_001.img - is the name of the file you want to "collect"
1s             - is the exposure time
200mm          - is the xtal-to-detector distance
param          - is one of: frames phi range exposure delta wedge dist energy wave
param                       flux mosaic beam_size xtal_size
value          - is the desired value for param

example:
mlfsom.com /data/${user}/test_1_001.img frames=100 1s 12657eV 200mm

EOF
exit 9
endif





# get crystal parameters from input file
echo "head" | mtzdump hklin $mtzfile >! ${tempfile}mtzdump
set RESO  = `awk '/Resolution Range/{getline;getline;print $6}' ${tempfile}mtzdump`
if($?USER_res) then
    set RESO = "$USER_res"
endif
set CELL  = `awk '/Cell Dimensions/{getline;getline;print}' ${tempfile}mtzdump`
set SGnum = `awk '/Space group/{print $NF+0}' ${tempfile}mtzdump`
set ASU_per_CELL = `awk -v SGnum=$SGnum '$1 == SGnum {print $2}' $CLIBD/symop.lib |& head -1`


# get all the names of "F" s
cat ${tempfile}mtzdump |\
awk '/ H K L /{for(i=1;i<=NF;++i)c[i]=$i} \
     / H H H /{for(i=1;i<=NF;++i)if($i=="G") print c[i]}' |\
cat >! ${tempfile}Gs
cat ${tempfile}mtzdump |\
awk '/ H K L /{for(i=1;i<=NF;++i)c[i]=$i} \
     / H H H /{for(i=1;i<=NF;++i)if($i=="F") print c[i]}' |\
cat >! ${tempfile}Fs
cat ${tempfile}mtzdump |\
awk '/ H K L /{for(i=1;i<=NF;++i)c[i]=$i} \
     / H H H /{for(i=1;i<=NF;++i)if($i=="D") print c[i]}' |\
cat >! ${tempfile}DANOs

# have a look at the data columns
set Fmean
set DANO = `head -1 ${tempfile}DANOs`

set Fplus = `awk 'toupper($0)=="FPLUS" || toupper($0)=="F(+)" || toupper($0)=="F+"' ${tempfile}Gs | head -1`
set Fminus = `awk 'toupper($0)=="FMINUS" || toupper($0)=="F(-)" || toupper($0)=="F-"' ${tempfile}Gs | head -1`
if("$Fplus" == "" || "Fminus" == "") then
    set Fplus = `awk 'toupper($0)=="FPLUS" || toupper($0)=="F(+)" || toupper($0)=="F+"' ${tempfile}Fs | head -1`
    set Fminus = `awk 'toupper($0)=="FMINUS" || toupper($0)=="F(-)" || toupper($0)=="F-"' ${tempfile}Fs | head -1`
endif
if("$Fplus" == "" || "Fminus" == "") then
    set Fplus = ""
    set Fmean = `awk 'toupper($0)=="FMEAN"' ${tempfile}Fs | head -1`
endif
if("$Fmean" == "") then
    set Fmean = `head -1 ${tempfile}Fs`
endif
if("$Fmean" == "" && "$Fplus" == "") then
    set BAD = "no Fs in $mtzfile"
    goto exit
endif

matthews_coef << end_mat >! ${tempfile}.log
CELL $CELL
SYMM P1
MOLW 1
NMOL 1
END
end_mat
set cellvolume = `awk '/Cell volume:/{print $NF}' ${tempfile}.log`

# calculate volume of xtal illuminated by the beam
set xtal_volume = `echo $xtal_size $beam_size | awk '{x=y=z=$1} $1>$2{x=y=$2} END{print x*y*z}'`

# number of unit cells in the beam
set unit_cells = `echo $xtal_volume $cellvolume | awk '{print ($1*1e-18)/($2*1e-30)}'`
set ASUs = `echo $unit_cells $ASU_per_CELL | awk '{print $1*$2}'`

# calculate the scale factor from Darwin's formula (relating I to |F|^2)
echo $brightness $exposure $osc $xtal_volume $cellvolume $wavelength |\
 awk '{PI=atan2(1,0)*2;dtr=PI/180;e=1.60217733e-19;m=9.1093897e-31;\
 c=299792458;h=6.6260755e-34;epsilon0=8.854187e-12;\
 re=(e^2)/(4*PI*epsilon0*m*c^2);\
 brightness=$1/(1e-6)^2;exposure=$2;osc=$3*dtr;xtalvolume=($4)*1e-18;cellvolume=$5*1e-30;lambda=$6*1e-10;\
 photonenergy=(h*c)/lambda;omega=osc/exposure;\
 powerdensity=brightness*photonenergy;\
 print brightness*exposure*(re)^2;\
 print brightness*(lambda^3)/omega*(re)^2*xtalvolume/(cellvolume^2)}' |\
cat >! ${tempfile}darwin
set darwin = `tail -1 ${tempfile}darwin`
set thomson = `tail -2 ${tempfile}darwin | head -1`

# full spot volume     (photons)     = darwin * lp * |F|^2
# incoherent intensity (photons/m^2) = thomson * electrons * (1+cos^2)/2




# make a bunch of Gaussian deviates
echo $frames $osc $mosaic |\
awk 'BEGIN{srand()} {for(i=1;i<=50*$1*(($2+1)/($3+1));++i){\
  rsq=0;\
  while((rsq>=1)||(rsq==0)){\
    x=2*rand()-1;y=2*rand()-1;rsq=x*x+y*y;\
  }fac=sqrt(-2*log(rsq)/rsq);\
  print x*fac}}' |\
cat >! ${tempfile}.gaussian_deviates
set gaussdev_index = 0



#################################################################################
# report

cat << EOF
Fs     $mtzfile

files  ${outprefix}_###.img
start  $firstbatch
end    $lastnum

cell   $CELL
mosaic $mosaic
rock   $rocking_curve

phi0   $phi0
phiend $phiend
range  $range
osc    $osc
frames $frames
nwaves $nwaves

dist   $distance
wave   $wavelength
misset $missets
matrix $Amatrix
beam   $beam_center
spindle_swing   $spindle_swing

xtal   $xtal_size um
EOF

set nbg = 0
foreach material ( $nonxtal_materials )
    @ nbg = ( $nbg + 1 )
    set thick = `echo $nonxtal_thickness[$nbg] | awk '$1+0>1e-3{print $1+0}'`
    if("$thick" != "") then
	echo "$material $thick" | awk '{printf "%-7s %s um\n",$1,$2}'
    endif
end

if("$icering_intensity" != "0") then
    echo "ice    ${icering_intensity}photons ${icering_width}pix"
endif
if(-e "$decayedmtz") then
    echo "half-dose $decay_halfdose MGy"
    echo "dose ratio: $decay_doseratio ph/um^2/Gy"
endif

echo ""

predict:
test -e preds.txt
if("$status" == "0") then
    echo "using ./preds.txt"
    rm -f  ${tempfile}_preds.txt >& /dev/null
    ln -sf `pwd`/preds.txt ${tempfile}_preds.txt
    goto expand
endif
if(-e ${tempfile}_preds.txt && $?CACHE) goto expand
###############################################################################
# record size for fake images
set record = `echo $dimx $bitsperpixel | awk '{print $1*$2/8}'`


# calculate size of detector area
set width = `echo $dimx $pixel | awk '{print $1*$2}'`

which hkl2XYphi.awk
echo "generating all spot predictions"
cat << EOF | tee ${tempfile}predin.txt
RESOLUTION    $RESO
WAVELENGTH    $wavelength
DISPERSION    $dispersion
DIVERGENCE    $divergence
MOSAIC        $mosaic
SPEEDUP_PIX   $speedup_pix
SPEEDUP_PHI   $speedup_phi
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
echo "monitoring a few: "
echo "   h    k    l        phi     delta-phi         L*p/A/y   subspots     Xdet        Ydet           tilt        width1      width2"
./hkl2XYphi.awk ${tempfile}predin.txt | tee ${tempfile}_preds.txt | awk '$1!=last{print;last=$1}'
# format: h k l   phi delta_phi   LP subspots  Xpix Ypix   psi Rsize Tsize
set preds = `cat ${tempfile}_preds.txt | wc -l`
echo ""
msdate.com $start_time
echo "s to complete $preds predictions"
set start_time = `msdate.com | awk '{print $NF}'`



#################################################################################
# now we have a look at the Fs

expand:
if(-e ${tempfile}pristine_Fsqr.hkl && $?CACHE) goto vignettemask

if(-e ${tempfile}_preds_1.XYI) then
    # this cache is no longer valid
    rm -f ${tempfile}_preds_*.XYI >& /dev/null
endif

if("$Fplus" == "" && "$Fmean" != "" && "$DANO" == "") then
 # no anomalous difference
 echo "getting $Fmean as Fplus and Fminus out of $mtzfile"
 cad hklin1 $mtzfile hklin2 $mtzfile \
     hklout ${tempfile}Fplus_Fminus.mtz <<EOF > /dev/null
labin file 1 E1=$Fmean
labou file 1 E1=Fplus
labin file 2 E1=$Fmean
labou file 2 E1=Fminus
EOF
endif

if("$Fplus" == "" && "$Fmean" != "" && "$DANO" != "") then
 # dump Fplus and Fminus
 echo "getting $Fmean and $DANO out of $mtzfile"
 sftools << EOF > /dev/null
read $mtzfile
calc COL Fplus = COL $Fmean COL $DANO 2 / +
calc COL Fminus = COL $Fmean COL $DANO 2 / -
write ${tempfile}Fplus_Fminus.mtz COL Fplus Fminus
y
stop
y
EOF

  if(-e "$decayedmtz") then
  # dump Fplus and Fminus
  echo "getting $Fmean and $DANO out of $decayedmtz"
  sftools << EOF > /dev/null
read $decayedmtz
calc COL Fplus = COL $Fmean COL $DANO 2 / +
calc COL Fminus = COL $Fmean COL $DANO 2 / -
write ${tempfile}decayed_Fplus_Fminus.mtz COL Fplus Fminus
y
stop
y
EOF
  endif
endif

if("$Fplus" != "" && "$Fminus" != "") then
echo "getting $Fplus and $Fminus out of $mtzfile"
cad hklin1 $mtzfile hklout ${tempfile}Fplus_Fminus.mtz <<EOF > /dev/null
labin file 1 E1=$Fplus E2=$Fminus
labou file 1 E1=Fplus  E2=Fminus
EOF

  if(-e "$decayedmtz") then
    echo "getting $Fplus and $Fminus out of $decayedmtz"
cad hklin1 $decayedmtz hklout ${tempfile}decayed_Fplus_Fminus.mtz <<EOF > /dev/null
labin file 1 E1=$Fplus E2=$Fminus
labou file 1 E1=Fplus  E2=Fminus
EOF

  endif

endif

#mtzutils hklin1 $mtzfile hklout ${tempfile}Fplus_Fminus.mtz << EOF > /dev/null
#include Fplus Fminus
#EOF

if(! -e ${tempfile}Fplus_Fminus.mtz) then
    set BAD = "unable to extract Fs"
    goto exit
endif


# figure out which hkls are F+ or F-
unique hklout ${tempfile}unique.mtz << EOF >&! ${tempfile}unique.log
RESO 1000
SYMM $SGnum
CELL $CELL
EOF
set Fplus_ops = `awk '/Bijvoet positive/{getline;while(NF){getline;for(i=2;i<NF;i+=2)print $(i+1)}}' ${tempfile}unique.log`
set Fminus_ops = `awk '/Bijvoet negative/{getline;while(NF){getline;for(i=2;i<NF;i+=2)print $(i+1)}}' ${tempfile}unique.log`



# loop over symmetry operators to get full reciprocal-space unit cell
echo -n "" >! ${tempfile}pristine_Fsqr.hkl
echo -n "" >! ${tempfile}decayed_Fsqr.hkl
foreach op ( $Fplus_ops )

reindex hklin ${tempfile}Fplus_Fminus.mtz hklout ${tempfile}asu.mtz << EOF > /dev/null
noreduce
reindex HKL $op
EOF

echo "extracting F+ and F- for $op"
awk 'BEGIN{print "lreso\nnref -1"}' |\
 mtzdump hklin ${tempfile}asu.mtz |\
awk '/ H K L /{for(i=1;i<=NF;++i){if($i=="Fplus")iFplus=i;if($i=="Fminus")iFminus=i}}\
  /LIST OF REFLECTIONS/{p=1} /SUMMARY/{p=0}\
  p && NF==6{h=$1;k=$2;l=$3;stol=sqrt($4/4);Fplus=$(iFplus+1);Fminus=$(iFminus+1);\
  print "HKL",h,k,l,stol,Fplus*Fplus;\
  print "HKL",-h,-k,-l,stol,Fminus*Fminus}' |\
awk '$6>0' >> ${tempfile}pristine_Fsqr.hkl
# format "HKL"  H  K  L sin(theta)/lambda |F+|^2
# format "HKL" -H -K -L sin(theta)/lambda |F-|^2


if(-e "$decayedmtz") then
reindex hklin ${tempfile}decayed_Fplus_Fminus.mtz hklout ${tempfile}asu.mtz << EOF > /dev/null
noreduce
reindex HKL $op
EOF

echo "extracting F+ and F- for $op from $decayedmtz"
awk 'BEGIN{print "lreso\nnref -1"}' |\
 mtzdump hklin ${tempfile}asu.mtz |\
awk '/ H K L /{for(i=1;i<=NF;++i){if($i=="Fplus")iFplus=i;if($i=="Fminus")iFminus=i}}\
  /LIST OF REFLECTIONS/{p=1} /SUMMARY/{p=0}\
  p && NF==6{h=$1;k=$2;l=$3;stol=sqrt($4/4);Fplus=$(iFplus+1);Fminus=$(iFminus+1);\
  print "HKL",h,k,l,stol,Fplus*Fplus;\
  print "HKL",-h,-k,-l,stol,Fminus*Fminus}' |\
awk '$6>0' >> ${tempfile}decayed_Fsqr.hkl
# format "HKL"  H  K  L sin(theta)/lambda |F+|^2
# format "HKL" -H -K -L sin(theta)/lambda |F-|^2
endif

end

msdate.com $start_time
echo "s to complete expansion"
set start_time = `msdate.com | awk '{print $NF}'`


vignettemask:
#########################################################################################

if(-e "${tempfile}_vignette.f2d" && $?CACHE) goto ripplemask
echo "creating vignette mask..."
# vignette effect for each module
# FIT2D's gaussian is exp(-( x^2/(2*xsig^2) + y^2/(2*ysig^2) ))
# so input "sigma" of fwhm/(2*sqrt(2*log(2)))
set vignette_sigma = `echo $vignette $dimx $windowpanes_x | awk '$3<=0{$3=1} $1+0<=0{print 999999;exit} {print ($2/$3)*sqrt(2)/sqrt(-log(1-$1)/log(2)) /(2*sqrt(2*log(2)))}'`
if("$vignette_sigma" == "0") set vignette_sigma = 999999;
#echo "vignette_sigma = $vignette_sigma pixels"

cat << EOF >! ${tempfile}.vignette
$dimx
$dimy
N
CREAT
$dimx
$dimy
EOF


foreach module ( `echo $dimx $windowpanes_x $dimy $windowpanes_y | awk '{xmodule=int($1/$2);ymodule=int($3/$4);for(x=0;x<$1;x+=xmodule)for(y=0;y<$1;y+=ymodule) if(x+xmodule<=$1&&y+ymodule<=$3) print x"_"y"_"xmodule"_"ymodule}'` )
    set x    = `echo $module | awk -F "_" '{print $1}'`
    set y    = `echo $module | awk -F "_" '{print $2}'`
    set xedge = `echo $module | awk -F "_" '{print $3}'`
    set yedge = `echo $module | awk -F "_" '{print $4}'`
    set ROI_x_lower = `echo $x 1     | awk '{print $1+$2}'`
    set ROI_y_lower = `echo $y 1     | awk '{print $1+$2}'`
    set ROI_x_upper = `echo $x $xedge | awk '{print $1+$2}'`
    set ROI_y_upper = `echo $y $yedge | awk '{print $1+$2}'`
    set center_x = `echo $x $xedge | awk '{print $1+$2/2}'`
    set center_y = `echo $y $yedge | awk '{print $1+$2/2}'`
cat << EOF >> ${tempfile}.vignette
ROI
$ROI_x_lower
$ROI_y_lower
$ROI_x_upper
$ROI_y_upper
GAUSS
$center_x
$center_y
1.0
0
$vignette_sigma
$vignette_sigma
EOF
end
cat << EOF >> ${tempfile}.vignette
ROI
1
1
$dimx
$dimy
cadd
1e-3
cmult
$peak_dqe
stat
output
fit2d
${tempfile}_vignette.f2d
EXIT
Y
EOF

cat ${tempfile}.vignette |\
fit2d -nographics >! ${tempfile}.vignette.log


ripplemask:
#########################################################################################

if(-e "${tempfile}_ripplemask.f2d" && $?CACHE) goto batchloop

if(-e ripplemask.f2d) then
    echo "using ./ripplemask.f2d"
    ln -sf `pwd`/ripplemask.f2d ${tempfile}_ripplemask.f2d
    goto batchloop
endif
echo "creating ${ripplenoise}% ripple mask with ~${ripplenoise_size}-pixel ripples"

# use the same "ripple mask" for all images
set ripplenoise_fac = `echo $ripplenoise | awk '{print $1/100}'`

cat << EOF >! ${tempfile}.ripplemask
$dimx
$dimy
N
CREAT
$dimx
$dimy
cadd
10000
poisson
cadd
-10000
blur
$ripplenoise_size
$ripplenoise_size
exchange
blur
$ripplenoise_size
$ripplenoise_size
exchange
blur
$ripplenoise_size
$ripplenoise_size
exchange
blur
$ripplenoise_size
$ripplenoise_size
exchange
stat
cdiv
##SIGMA
cmult
$ripplenoise_fac
cadd
1
stat
output
fit2d
${tempfile}_ripplemask.f2d
EXIT
Y
EOF

cat ${tempfile}.ripplemask |\
fit2d -nographics >! ${tempfile}.ripplemask.log


batchloop:
set >! ${tempfile}settings.txt
#################################################################################
foreach frame ( `seq 1 1 $frames` )

set batch = `echo $firstbatch $frame | awk '{print $1+$2-1}'`
set num = `echo $batch | awk '{printf "%03d",$1}'`
set phi = `echo "$phi0 $osc $batch $firstbatch" | awk '{print $1+$2*($3-$4)}'`
set phiend = `echo "$phi0 $osc $batch $firstbatch" | awk '{print $1+$2*($3-$4+1)}'`


set outfile = `echo $outprefix $batch | awk '{printf "%s_%03d.img",$1,$2}'`
if(-e $outfile && $?CACHE) then
    echo "$outfile" already exists
    if("$CACHE" == "images") continue
endif
# make sure no other jobs try to make this
if($?CACHE) then
    if("$CACHE" == "images") then
       echo "$HOST $$" >! $outfile
    endif
endif


# skip every so many?
set test = `echo $frame $skip_step $skip_start | awk '{print ($1*$2*$3==0 || $1%$2==($3-1))}'`
if(! $test) then
    echo "skipping $outfile"
    continue
endif


echo "batch $batch	phi: $phi - $phiend"

# allow for frame-by-frame control of various things
if(-e flux_schedule.txt) then
    # format: batch flux xtalsize
    set line = `awk -v batch=$batch '$1+0==batch+0' flux_schedule.txt`
    set test = `echo $line | awk '$2+0>0{print $2+0}'`
    if("$test" != "") then
	set flux = $test
	set test = `echo $flux | awk '{printf "%.5g",$1}'`
	echo "new flux: $test ph/s"
	set brightness = `echo $flux $beam_size | awk '{print $1/$2/$2}'`
    endif
    set test = `echo $line | awk '$3+0>0{print $3+0}'`
    if("$test" != "") then
	set xtal_size = "$test"
	echo "new xtal size: $xtal_size um"

    endif
endif


# calculate progression of radiation damage
set total_exposure_time = `echo $batch $exposure | awk '{print $1*$2}'`
# compute the dose in MGray
set dose = `echo $total_exposure_time $brightness $decay_doseratio $starting_dose | awk '{print $1*$2/$3/1e6 + $4}'`
# compute the fractional progression of an exponential decay
set frac_decayed = `echo $dose $decay_halfdose | awk '{print 1-exp(-log(2)*$1/$2)}'`
# s * ph/s/um^2 / ph/um^2/Gy / Gy = fraction
# compute additional B factor from dose
set decay_A = `echo $dose $decay_AperMGray | awk '{print $1*$2}'`

echo "total exposure is now $total_exposure_time s"
echo "brightness is now $brightness ph/s/um^2"
echo "total dose is now $dose MGy"
echo "damage A factor is now $decay_A"
echo "damage reaction is now $frac_decayed fractional complete"

# SHUTTER JITTER
# calculate deviations to phi range given that shutter is not perfect
@ gaussdev_index = ( $gaussdev_index + 1)
set gaussdev = `head -$gaussdev_index ${tempfile}.gaussian_deviates | tail -1`
set noisy_phi = `echo $phi $osc $exposure $shutter_jitter $gaussdev | awk '{print $1+ $2/$3*$4*1e-3*$5}'`
@ gaussdev_index = ( $gaussdev_index + 1)
set gaussdev = `head -$gaussdev_index ${tempfile}.gaussian_deviates | tail -1`
set noisy_phiend = `echo $phiend $osc $exposure $shutter_jitter $shutter_gap $gaussdev | awk '{print $1+ $2/$3*(($4*$6)-$5)*1e-3}'`
# assume motor speed is stable?
set angular_speed  = `echo $osc $exposure | awk '{print $1/$2}'`
set noisy_osc      = `echo $noisy_phiend $noisy_phi | awk '{print ($1-$2)}'`
set noisy_exposure = `echo $noisy_osc $angular_speed | awk '{print $1/$2}'`

# BEAM FLICKER
# calculate "effecetive" brightness given that the beam is noisy...
@ gaussdev_index = ( $gaussdev_index + 1)
set gaussdev = `head -$gaussdev_index ${tempfile}.gaussian_deviates | tail -1`
set noisy_brightness = `echo $brightness $beam_flicker $gaussdev $noisy_exposure | awk '{print $1+$1*($2/100)*$3/sqrt($4)}'`


# forget all this noise if user does not want it
if($?NO_NOISE) then
    set noisy_brightness = $brightness
    set noisy_exposure = $exposure
    set noisy_phi      = $phi
    set noisy_phiend   = $phiend
    set noisy_osc      = $osc
    set subranges      = 1
    set beam_flicker   = 0
endif
#echo "noisy brightness: $noisy_brightness ph/s/um^2"

# calculate volume of xtal illuminated by the beam (in um^3)
set xtal_volume = `echo $xtal_size $beam_size | awk '{x=y=z=$1} $1>$2{x=y=$2} END{print x*y*z}'`

# number of unit cells in the beam
set unit_cells = `echo $xtal_volume $cellvolume | awk '{print ($1*1e-18)/($2*1e-30)}'`
set ASUs = `echo $unit_cells $ASU_per_CELL | awk '{print $1*$2}'`

# calculate the scale factor from Darwin's formula (relating I to |F|^2)
# noisy_brightness will be implemented in the partiality calculation
echo $brightness $noisy_exposure $noisy_osc $xtal_volume $cellvolume $wavelength |\
 awk '{PI=atan2(1,0)*2;dtr=PI/180;e=1.60217733e-19;m=9.1093897e-31;\
 c=299792458;h=6.6260755e-34;epsilon0=8.854187e-12;\
 re=(e^2)/(4*PI*epsilon0*m*c^2);\
 brightness=$1/(1e-6)^2;exposure=$2;osc=$3*dtr;xtalvolume=($4)*1e-18;cellvolume=$5*1e-30;lambda=$6*1e-10;\
 photonenergy=(h*c)/lambda;omega=osc/exposure;\
 powerdensity=brightness*photonenergy;\
 print brightness*exposure*(re)^2;\
 print brightness*(lambda^3)/omega*(re)^2*xtalvolume/(cellvolume^2)}' |\
cat >! ${tempfile}darwin
set darwin = `tail -1 ${tempfile}darwin`
set thomson = `tail -2 ${tempfile}darwin | head -1`

# full spot volume     (photons)     = darwin * lp * |F|^2
# incoherent intensity (photons/m^2) = thomson * electrons * (1+cos^2)/2


# define the integral of the mosaic spread (rocking curve)
cat << EOF >! ${tempfile}_rockfunc
s=sqrt(8.0*log(2))
gauss_rock(x) = norm(x*s)
lorentz_rock(x) = atan(x*2)/pi
fwhm=1.76274717126295
tanh_rock(x) = tanh(x*fwhm)/2
square_rock(x) = ((x<-0.5)?0:((x>0.5)?1:(x+0.5)))
sr = 10
slope=sqrt(4*log(2)/pi)
area=1.0*(1-1.0/sr)+1.0/sr/slope
area=1-1.0/sr+1.0/sr/slope
tophat_rock(x) = ((x<-(0.5-0.5/sr))?gauss_rock((x+0.5-0.5/sr)*sr)/sr/slope:((x>0.5-0.5/sr)?(gauss_rock((x-0.5+0.5/sr)*sr)/sr/slope+(1.-1./sr)):(0.5/sr/slope+(x+(0.5-0.5/sr)))))/(1-1./sr+1./sr/slope)
disk_rock(x) = (x<=-0.5?0:(x>=0.5?1:(2*x*sqrt(1-4*x**2)+asin(2*x))/pi+0.5))
EOF
if("$rocking_curve" == "tanh") then
    set etamax = 7
    echo "rock(x) = tanh_rock(x)" >> ${tempfile}_rockfunc
endif
if("$rocking_curve" == "gauss") then
    set etamax = 3
    echo "rock(x) = gauss_rock(x)" >> ${tempfile}_rockfunc
endif
if("$rocking_curve" == "lorentz") then
    set etamax = 1e9
    echo "rock(x) = lorentz_rock(x)" >> ${tempfile}_rockfunc
endif
if("$rocking_curve" == "square") then
    set etamax = 1.1
    echo "rock(x) = square_rock(x)" >> ${tempfile}_rockfunc
endif
if("$rocking_curve" == "tophat") then
    set etamax = 2
    echo "rock(x) = tophat_rock(x)" >> ${tempfile}_rockfunc
endif
if("$rocking_curve" =~ dis[ck]) then
    set etamax = 1.1
    echo "rock(x) = disk_rock(x)" >> ${tempfile}_rockfunc
endif


# break up phi range over several sub-ranges for computing flicker noise
set subranges = 10
if("$beam_flicker" == "0") set subranges = 1
echo "using $subranges phi subranges"
cat << EOF >> ${tempfile}_rockfunc
part(start,end)=(\
EOF
cat ${tempfile}.gaussian_deviates |\
awk -v i=$gaussdev_index -v subranges=$subranges 'NR>i{print;++n} n>=subranges{exit}' |\
awk -v subranges=$subranges -v flicker0=$beam_flicker -v time=$noisy_exposure '\
   {gaussdev=$1;print (1 + flicker0/100*gaussdev/sqrt(time/subranges))}' |\
cat >! ${tempfile}.phiweights
@ gaussdev_index = ( $gaussdev_index + $subranges )

set overall_flicker = `awk '{++n;sum+=$1} END{print sum/n}' ${tempfile}.phiweights`
echo "overall flicker: $overall_flicker"
set noisy_darwin = `echo $overall_flicker $darwin | awk '{print $1*$2}'`
set noisy_thomson = `echo $overall_flicker $thomson | awk '{print $1*$2}'`

cat ${tempfile}.phiweights |\
awk -v subranges=$subranges '\
   {print "+"$1"*(rock(start+"n+0"./"subranges".*(end-start))-rock(start+"n+1"./"subranges".*(end-start)))\\";++n}' >> ${tempfile}_rockfunc
echo ")" >> ${tempfile}_rockfunc



if(-e ${tempfile}_preds_${batch}.XY && $?CACHE) goto FtoI
echo "extracting relavant predictions"


# compute partiality: using above integrated rocking curve
cat ${tempfile}_preds.txt |\
awk -v phistart=$noisy_phi -v phiend=$noisy_phiend \
    -v etamax=$etamax '\
 function degdiff(x,y) {return (x-y+36180)%360-180}\
 {h=$1;k=$2;l=$3;phi=$4;dphi=$5;lps=$6*$7;X=$8;Y=$9;\
  psi=$10;Rsize=$11;Tsize=$12}\
 dphi>0{etastart=degdiff(phistart,phi)/dphi;\
    etaend=degdiff(phiend,phi)/dphi;\
  if(etastart<-etamax)etastart=-etamax;\
  if(etastart> etamax)next;\
  if(etaend  <-etamax)next;\
  if(etaend  > etamax)etaend  = etamax;\
  if(etastart==etaend){next}\
  print X,Y,h,k,l,lps,etastart,etaend,phi,psi,Rsize,Tsize}' |\
tee ${tempfile}_etarange |\
awk '{print "print part("$8","$7")"}' |\
cat ${tempfile}_rockfunc - |\
gnuplot |&\
cat >! ${tempfile}_parts

# since gnuplot is slow with text I/O, re-assemble list here
cat ${tempfile}_parts ${tempfile}_etarange |\
awk 'NF==1{if($1!="0.0")part[NR]=$1+0;next}\
     {++n;if(part[n]+0>1e-10)print $1,$2,$3,$4,$5,$6,part[n],$9,$10,$11,$12}' |\
cat >! ${tempfile}_preds_${batch}.XY
# format: XDET YDET H K L Lorentz*polar*subspots/yield frac_recorded PHI psi Rsize Tsize


msdate.com $start_time
echo "s to extract `cat ${tempfile}_preds_${batch}.XY | wc -l` preds"
set start_time = `msdate.com | awk '{print $NF}'`

FtoI:
#if(-e ${tempfile}_preds_${batch}.XYI && $?CACHE) goto generate
if(-e "$decayedmtz") then
    echo "simulating decay fraction: $frac_decayed"
    # combine decayed data with undecayed data
    awk '{print $0,"DECAYED"}' ${tempfile}decayed_Fsqr.hkl |\
    cat - ${tempfile}pristine_Fsqr.hkl |\
    awk -v frac=$frac_decayed '/DECAYED/{Idecayed[$2,$3,$4]=$6;next} \
       /^HKL/{print $1,$2,$3,$4,$5,((1-frac)*sqrt($6)+frac*sqrt(Idecayed[$2,$3,$4]))^2}' |\
    cat >! ${tempfile}Fsqr.hkl

    msdate.com $start_time
    echo "s to simulate decay fraction $frac_decayed"
    set start_time = `msdate.com | awk '{print $NF}'`
else
    ln -sf ${tempfile}pristine_Fsqr.hkl ${tempfile}Fsqr.hkl
endif

set total_B = `echo $extra_B | awk '{print $1}'`
set total_A = `echo $extra_A $decay_A | awk '{print $1+$2}'`

# now we encode the Fsquared values into spot intensities
echo "calculating spot intensities"
cat ${tempfile}Fsqr.hkl ${tempfile}_preds_${batch}.XY |\
awk -v darwin=$darwin -v B=$total_B -v A=$total_A '\
 /^HKL/{stol=$5;e=2*A*stol+2*B*stol^2;if(e>100 || e<-100){next};\
  I[$2,$3,$4]=$6*exp(-e);next} \
  NF>=7 && $6*$7>0 {X=$1;Y=$2;h=$3;k=$4;l=$5;lp=$6;frac=$7;phi=$8;psi=$9;Rsize=$10;Tsize=$11;\
  if(I[h,k,l]=="") I[h,k,l]=I[-h,-k,-l];\
  if(I[h,k,l]=="") next;\
  print X,Y,darwin*I[h,k,l]/lp*frac,psi,Rsize,Tsize,h,k,l}' |\
cat >! ${tempfile}_preds_${batch}.XYI
# format: XDET YDET photons psi Rsize Tsize  h k l

msdate.com $start_time
echo "s to calculate `cat ${tempfile}_preds_${batch}.XYI | wc -l` intensities"
set start_time = `msdate.com | awk '{print $NF}'`

generate:
#################################################################################
# generate the images


echo "generating $outfile"


set SN = `echo $windowpanes_x $dimx $pixel | awk '{mods=$1*$1;w=$2*$3} mods==9{print "900";exit} w>200{print 440;exit} {print mods"00"}'`
cat << EOF >! ${tempfile}.header
{
HEADER_BYTES=  512;
DIM=2;
BYTE_ORDER=little_endian;
TYPE=unsigned_short;
SIZE1=$dimx;
SIZE2=$dimy;
BIN=2x2;
DETECTOR_SN=$SN;
BEAMLINE=FAKE;
PIXEL_SIZE=$pixel;
BEAM_CENTER_X=$adxv_center[1];
BEAM_CENTER_Y=$adxv_center[2];
DATE=`date`;
TIME=$exposure;
DISTANCE=$distance;
TWOTHETA=$twotheta;
OSC_RANGE=$osc;
PHI=$phi;
OSC_START=$phi;
AXIS=phi;
WAVELENGTH=$wavelength;
}
EOF

# opening FIT2D commands
cat << EOF >! ${tempfile}.start
$dimx
$dimy
N
CREAT
$dimx
$dimy
EOF

# closing FIT2D commands to write out binary integer image
set record = `echo $dimx $bitsperpixel | awk '{print $1*$2/8}'`
set maxvalue = `echo $bitsperpixel | awk '{print 2^$1-1}'`
cat << EOF >! ${tempfile}.end
OUTPUT
bin
${tempfile}fit2d$$.bin
${record}
1
N
0
$maxvalue
EXIT
Y
EOF


# add shot noise (after applying the vignette mask)
# add read noise (calculate as noise-equivalent power
# using (readnoise/EOgain)^2 to convert to "equivalent photon fog":
# this number of photons/pixel will deposit the same amount of noise as the read noise
set readnoise_photons = `echo $readnoise $EOgain | awk '{print ($1/$2)^2}'`
# this is the rms error (in photon equivalents/pixel) that the read noise will create
set readnoise_rmsphotons = `echo $readnoise $EOgain | awk '{print $1/$2}'`
set inv_gain = `echo $gain | awk '{print 1/$1}'`

# a finite point-spread will smooth out our pure Gaussian read noise model
# so, we need to crank up the noise here to combat the smoothing function
set readnoise_scale = `echo $psf $pixel | awk '{x=$1/$2;print (x/4.22914+1)^2.78836}'`

set readnoise_raw = `echo $readnoise_scale $readnoise_rmsphotons | awk '{print $1*$2}'`



cat << EOF >! ${tempfile}.noise
exchange
input
fit2d
${tempfile}_vignette.f2d
exchange
multiply
poisson
divide
exchange
creat
$dimx
$dimy
cadd
10000
poisson
cadd
-10000
cdiv
100
stat
cmult
$readnoise_raw
add
EOF


cat << EOF >! ${tempfile}.ripplenoise
exchange
input
fit2d
${tempfile}_ripplemask.f2d
multiply
stat
EOF


# point-spread function simulated by a Gussian image blur
if ("$psf" != "0") then
    # oversample the Gaussian by 10x since the shape of the PSF
    # probably varies significantly from one side of the pixel to the other
    echo $psf $pixel |\
    awk '{sigma=$1/$2;for(x=-1;x<=1;x+=0.1)for(y=-1;y<=1;y+=0.1){\
	printf "%.0f %.0f %s\n", x,y,exp(-log(2)*(x^2+y^2)/sigma^2)}}' |\
    awk '{printf "%.0f %.0f %s\n",$1+1e-3,$2+1e-3,$3}' |\
    awk '{++n;sum+=$3;coff[$1,$2]+=$3} END{norm=1/sum;\
	for(x=-1;x<=1;++x)for(y=-1;y<=1;++y) print coff[x,y]*norm}' |\
    awk 'BEGIN{\
	 print "spatial";\
	 print "y";\
	 print 4}\
        {print}\
     END{print "exchange"}' |\
    cat >! ${tempfile}.psf
else
   echo -n "" >! ${tempfile}.psf
endif
echo "stat" >> ${tempfile}.psf


# generate "window-pane" effect for tiled detectors
set seam = 0
echo -n "" >! ${tempfile}.windowpane
set modx = `echo $dimx $windowpanes_x $windowpane_width_x | awk '{print ($1-($2-1)*$3)/$2}'`
set mody = `echo $dimy $windowpanes_y $windowpane_width_y | awk '{print ($1-($2-1)*$3)/$2}'`
while ( $seam < $windowpanes_x )
@ seam = ( $seam + 1 )
set minpos = `echo $seam $modx $windowpane_width_x | awk '{print $1*$2+$3*($1-1)}'`
set maxpos = `echo $minpos $windowpane_width_x | awk '{print $1+$2-1}'`
set minpos = `echo $minpos $dimx | awk '$1>$2+0{$1=$2} {print $1+0}'`
set maxpos = `echo $maxpos $dimx | awk '$1>$2+0{$1=$2} {print $1+0}'`
if("$minpos" == "$maxpos") continue
cat << EOF >> ${tempfile}.windowpane
REGION
$minpos
0.01
$maxpos
$dimy
clear
cadd
$windowpane_value
EOF
end
set seam = 0
while ( $seam < $windowpanes_y )
@ seam = ( $seam + 1 )
set minpos = `echo $seam $mody $windowpane_width_y | awk '{print $1*$2+$3*($1-1)}'`
set maxpos = `echo $minpos $windowpane_width_y | awk '{print $1+$2-1}'`
set minpos = `echo $minpos $dimy | awk '$1>$2+0{$1=$2} {print $1+0}'`
set maxpos = `echo $maxpos $dimy | awk '$1>$2+0{$1=$2} {print $1+0}'`
if("$minpos" == "$maxpos") continue
cat << EOF >> ${tempfile}.windowpane
REGION
0.01
$minpos
$dimx
$maxpos
clear
cadd
$windowpane_value
EOF
end
echo "full region" >> ${tempfile}.windowpane


# apply the final apparent ADU/photon gain
cat << EOF >! ${tempfile}.adc
CMULT
$gain
CADD
$adcoffset
EOF


#################################
# type of spot to plot

# FIT2D's gaussian is exp(-( x^2/(2*xsig^2) + y^2/(2*ysig^2) ))
# the volume under this curve is 2*pi
# so, to plot a Gaussian with a given volume and fwhm:
# give FIT2D a "width"  of fwhm/(2*sqrt(2*log(2)))
# give FIT2D a "height" of volume/(fwhm*fwhm)*4*log(2)/pi
# the volume of the gaussian will be: height*fwhm*fwhm*pi/(4*log(2))

set fit2d_spotwidth = `echo $spotwidth | awk '{print $1/(2*sqrt(2*log(2)))}'`

if( "$spot_type" == "simple_gauss" ) then

echo "plotting simple gaussian spots"
cat ${tempfile}_preds_${batch}.XYI |\
awk -v spotwidth=$spotwidth 'BEGIN{PI=2*atan2(1,0);\
  fws=1./(2*sqrt(2*log(2)));fhs=4*log(2)/PI}\
  NF>2 && $3>0{volume=$3;y=$1;x=$2;\
    height=volume/spotwidth/spotwidth*fhs;\
    print "GAUSS";print x; print y;\
    print height; print 0; print spotwidth*fws; print spotwidth*fws;\
}' >! ${tempfile}_${batch}.spots

endif

if( "$spot_type" == "gauss" ) then

echo "plotting gaussian spots"
set db_pix = `echo $beam_center $pixel | awk '{print $1/$3,$2/$3}'`
set dist_pix = `echo $distance $pixel | awk '{print $1/$2}'`
cat ${tempfile}_preds_${batch}.XYI |\
awk 'BEGIN{PI=2*atan2(1,0);fws=1./(2*sqrt(2*log(2)));fhs=4*log(2)/PI}\
  NF>2 && $3>0.1{y=$1;x=$2;volume=$3;\
    tilt=$4;radsize=$5;tansize=$6;\
    height=volume/radsize/tansize*fhs;\
    if(height<0.5)next;\
    print "GAUSS";print x; print y;\
    print height; print tilt; print radsize*fws; print tansize*fws;\
}' >! ${tempfile}_${batch}.spots

endif

if( "$spot_type" == "lg" ) then

# use a 2-gaussian aproximation to a 10:90 lorentzian-gaussian function
echo "plotting lorentzian-gaussian spots"
cat ${tempfile}_preds_${batch}.XYI |\
awk 'BEGIN{PI=2*atan2(1,0);fws=1./(2*sqrt(2*log(2)));fhs=4*log(2)/PI;\
  # fitted parameters to the lg function\
  h1=0.975343;h3=2.5e-4;h2=1-h1-h3;w1=0.987952;w2=2.92039;w3=10;\
  # calculate relative volumes for conservation below \
  v1=h1*w1*w1;v2=h2*w2*w2;v3=h3*w3*w3;v_total=v1+v2+v3;\
  frac1=v1/v_total;frac2=v2/v_total;frac3=v3/v_total;\
  }\
  NF>2 && $3>0.1{y=$1;x=$2;volume=$3;\
    tilt=$4;radsize=$5;tansize=$6;\
    # relative widths of gaussians \
    # larger components are less sensitive to spot shape \
    gw1r=radsize*w1;gw1t=tansize*w1;\
    gw2r=gw1r+w2;gw2t=gw1t+w2;\
    gw3r=gw1r+w3;gw3t=gw1t+w3;\
    # calculate peak heights that will conserve fractional volume of each component \
    gh1=volume*frac1/gw1r/gw1t*fhs;\
    gh2=volume*frac2/gw2r/gw2t*fhs;\
    gh3=volume*frac3/gw3r/gw3t*fhs;\
    print "GAUSS";print x; print y;\
    print gh1; print tilt; print gw1r*fws; print gw1t*fws;\
    print "GAUSS";print x; print y;\
    print gh2; print tilt; print gw2r*fws; print gw2t*fws;\
    print "GAUSS";print x; print y;\
    print gh3; print tilt; print gw3r*fws; print gw3t*fws;\
}' >! ${tempfile}_${batch}.spots

endif


if( "$spot_type" == "lg_vpsf" ) then

# use a 2-gaussian aproximation to a 10:90 lorentzian-gaussian function and a variable point spread
echo "plotting lorentzian-gaussian spots with variable PSF"
cat ${tempfile}_preds_${batch}.XYI |\
awk 'BEGIN{PI=2*atan2(1,0);f2d_sws=1/(2*sqrt(2*log(2)))\
  # fitted parameters to the lg function\
  h1=0.975343;h2=1-h1;w1=0.987952;w2=2.92039;vscale=0.861668934987687;\
  hscale=vscale/((f2d_sws)^2)/2/PI;\
  gw1=f2d_sws*w1;gw2=f2d_sws*w2}\
  NF>2 && $3>0.1{y=$1;x=$2;volume=$3;\
    tilt=$4;radsize=$5;tansize=$6;\
    modx=x%1024-512;mody=y%1024-512;rsqr=modx*modx+mody*mody;\
    blur=1+1*(2-2*exp(-rsqr/400000));\
    blurradsize=blur*radsize;\
    blurtansize=blur*tansize;\
    height=hscale*volume/blurradsize/blurtansize;\
    print "GAUSS";print x; print y;\
    print height*h1; print tilt; print blurradsize*gw1; print blurtansize*gw1;\
    print "GAUSS";print x; print y;\
    print height*h2; print tilt; print blurradsize*gw2; print blurtansize*gw2;\
}' >! ${tempfile}_${batch}.spots

endif


if( "$spot_type" == "lorentz" ) then

# use a 3-gaussian aproximation to lorentzian
echo "plotting lorentzian spots"
cat ${tempfile}_preds_${batch}.XYI |\
awk 'BEGIN{PI=2*atan2(1,0);f2d_sws=1/(2*sqrt(2*log(2)))}\
  NF>2 && $3>0.1{y=$1;x=$2;volume=$3;\
    tilt=$4;radsize=$5;tansize=$6;\
    height=volume/radsize/tansize;\
    print "GAUSS";print x; print y;\
    print height*0.52538; print tilt; print radsize*f2d_sws*0.68902; print tansize*f2d_sws*0.68902;\
    print "GAUSS";print x; print y;\
    print height*0.388629; print tilt; print radsize*f2d_sws*1.60594; print tansize*f2d_sws*1.60594;\
    print "GAUSS";print x; print y;\
    print height*0.0788759; print tilt; print radsize*f2d_sws*5.09314; print tansize*f2d_sws*5.09314;\
}' >! ${tempfile}_${batch}.spots

endif

if( "$spot_type" == "pinpoints" || "$spot_type" == "square" ) then

cat ${tempfile}_preds_${batch}.XYI |\
awk -v spotwidth=$spotwidth 'BEGIN{spotwidth-=1} NF>2 && $3>0{\
  volume=$3;y=$1;x=$2;\
  print int(x-0.5*spotwidth),int(y-0.5*spotwidth),int(x+0.5*spotwidth),int(y+0.5*spotwidth),volume}' |\
awk -v dimx=$dimx -v dimy=$dimy '\
  $1>0 && $2>0 && $3<dimx && $4<dimy{volume=$3;y=$1;x=$2;\
    print "REGION";print $1; print $2; print $3; print $4;\
    print "CADD"; print int($5/(1+$3-$1)/(1+$4-$2)); print "full region";\
}' >! ${tempfile}_${batch}.spots

endif


if( "$spot_type" == "fiber" ) then

# use an n-gaussian aproximation to fiberoptic psf(x,y,g) = g/(2*pi)*sqrt(g^2+x^2+y^2)^-3
echo "plotting spots with fiberoptic PSF g_center=$fiber_spread g_corner/g_center=$fiber_length_ratio"
cat ${tempfile}_preds_${batch}.XYI |\
awk -v pixel=$pixel -v g0=$fiber_spread -v gratio=$fiber_length_ratio \
  'BEGIN{PI=2*atan2(1,0);fws=1./(2*sqrt(2*log(2)));fhs=4*log(2)/PI;\
         minpix=0.1;g0/=1000}\
  NF>2 && $3>0.1{y=$1;x=$2;photons=$3;\
    tilt=$4;radsize=$5;tansize=$6;\
    # compute position on fiber module \
    modx=x%1024-512;mody=y%1024-512;rsqr=modx*modx+mody*mody;\
    # make g increase to g*gratio at corners \
    g=g0*(1-(1-gratio)*(2-2*exp(-rsqr/760000)));\
    # calculate needed gaussians \
    # least-squares fit to PSF \
    h[1] = 0.401006; w[1] = 1.42454;\
    h[2] = 0.0870611; w[2] = 3.5227;\
    h[3] = 0.00662917; w[3] = 8.89032;\
    h[4] = 0.000422407; w[4] = 22.5083;\
    h[5] = 2.61459e-05; w[5] = 57.0133;\
    h[6] = 1.61026e-06; w[6] = 144.431;\
    h[7] = 9.90565e-08; w[7] = 365.925;\
    h[8] = 6.09001e-09; w[8] = 927.213;\
    h[9] = 3.74256e-10; w[9] = 2349.73;\
    h[10] = 2.29924e-11; w[10] = 5955.03;\
    h[11] = 1.41247e-12; w[11] = 15092.5;\
    h[12] = 8.67552e-14; w[12] = 38251.7;\
    h[13] = 5.32877e-15; w[13] = 96925.9;\
    h[14] = 3.27697e-16; w[14] = 244805.0;\
    h[15] = 2.24084e-17; w[15] = 631050.0;\
    # make sure we accomodate REALLY bright spots \
    for(i=15;photons*h[i-1]/PI*(g/pixel)^2>minpix;++i){\
     w[i]=2.53219*w[i-1];h[i]=4.74216*(w[i])^-3;\
    };\
    # see how many terms we need \
    for(i=1;photons*h[i]/PI*(g/pixel)^2>minpix;++i){};\
    n=i;if(n<1)n=1;\
    # calculate height and volume fraction of each component \
    h0=cf=0;\
    for(i=1;i<=n;++i){\
      h0+=h[i];frac[i]=h[i]*w[i]*w[i]/log(16);cf+=frac[i];cum_frac[i]=cf;\
    };\
    # compensate for correct total height \
    extrah=0.5-h0;\
    for(i=1;i<=n;++i){h[i]+=extrah*h[i]/h0};\
    # re-calculate height and volume fraction of each component \
    h0=cf=0;\
    for(i=1;i<=n;++i){\
      h0+=h[i];frac[i]=h[i]*w[i]*w[i]/log(16);cf+=frac[i];cum_frac[i]=cf;\
    };\
    # compensate lowest Gaussian for unit total volume \
    wnsqr=(1-cum_frac[n-1])*log(16)/h[n];\
    hn=(1-cum_frac[n-1])*log(16)/w[n]/w[n];\
    if(wnsqr>0){w[n]=sqrt(wnsqr)}else{h[n]=hn};\
    frac[n]=h[n]*w[n]*w[n]/log(16);\
    # stretch to given value of g and convert to pixel units \
    for(i=1;i<=n;++i){\
	h[i]=h[i]/PI/(g/pixel)^2;w[i]=w[i]*g/pixel;\
    }\
    for(i=1;i<=n;++i){\
      # convolute the spot with each component of the PSF \
      gwr[i]=sqrt(radsize^2+w[i]^2);gwt[i]=sqrt(tansize^2+w[i]^2);\
      # calculate peak heights that will conserve fractional volume of each component \
      gh[i]=photons*frac[i]/gwr[i]/gwt[i];\
      print "GAUSS";print x; print y;\
      print gh[i]*fhs; print tilt; print gwr[i]*fws; print gwt[i]*fws;\
    }\
}' >! ${tempfile}_${batch}.spots

endif




##############################################

# make some SAXS scatter?
echo "250 $SAXS_intensity $wavelength $distance $beam_center $pixel $noisy_exposure" |\
awk 'BEGIN{PI=2*atan2(1,0)} {\
   d=$1;height=$2;lambda=$3;dist=$4;x=$6/$7;y=$5/$7;pix=$7;exposure=$8;\
   sintheta=lambda/2/d;theta=atan2(sintheta,sqrt(1-sintheta^2));\
   r=dist*sin(2*theta)/pix;\
   print "GAUSS"; print x; print y;\
   print height*exposure; print 0; print r; print r;}' |\
cat >! ${tempfile}.saxs


# make a beamstop shadow
echo "$beamstop_size $beamstop_dist $distance $beam_center $pixel" |\
awk 'BEGIN{srand();f2d_sws=1/(2*sqrt(2*log(2)))} {\
   bs_size=$1;bs_dist=$2;detector_dist=$3;x=$5/$6;y=$4/$6;pix=$6;\
   shadow_radius=bs_size/bs_dist*detector_dist/pix*f2d_sws;\
   x+=2*rand()-1;y+=2*rand()-1;\
   print "STORE";\
   print "clear";\
   print "GAUSS"; print x; print y;\
   print 1; print 0; print shadow_radius; print 0.9*shadow_radius;\
   print "GAUSS"; print x; print y;\
   print 10000; print 0; print shadow_radius/4; print 0.9*shadow_radius/4;\
   print "CADD";\
   print "0.9";\
   print "THRESHOLD";\
   print "1";\
   print "9900";\
   print "EXCHANGE";\
   print "divide";\
}' |\
cat >! ${tempfile}.beamstop


# make the direct-beam spot
echo "1000 $divergence $distance  $beam_center  $pixel" |\
awk '{height=$1;Hdiv=$2/57.3;Vdiv=$3/57.3;dist=$4;pixel=$7;y=$5/pixel;x=$6/pixel;\
    f2d_sws=1/(2*sqrt(2*log(2)))\
    hsize=Hdiv*dist/pixel;vsize=Vdiv*dist/pixel;\
    for(hstep=-hsize/2;hstep<=hsize/2;++hstep){\
    for(vstep=-vsize/2;vstep<=vsize/2;++vstep){\
        print "GAUSS";print x+hstep; print y+vstep;\
        print height; print 0; print 1*f2d_sws; print 1*f2d_sws;\
    }}\
}' >! ${tempfile}.directbeam











# incoherent scattering
# elastic (diffuse) scatter from xtal
# AND Compton form xtal
# fluorescence from xtal?


# we need to know the B factors from the PDB to get the diffuse scattering strength
# we also need the atomic composition of the unit cell to get accurate Compton
cat $pdbfile |\
awk '/^ATOM/ || /^HETAT/{Ee=substr($0, 13, 2);B=substr($0, 61, 6)+0;\
  print Ee,B}' |\
awk -v total_A=$total_A -v total_B=$total_B '{Ee=$1;B=$2+total_B;\
 ++count[Ee];\
 for(stol=0;stol<=2.0;stol+=0.01){\
    exponent=2*A*stol+2*B*stol^2;\
    if(exponent>100){fac=0}else{fac=exp(-exponent)};\
    diff_scat_atoms[Ee" "stol]+=(1-fac);\
  }\
} END{\
    for(Ee_stol in diff_scat_atoms) {\
	split(Ee_stol,w);Ee=w[1];\
	print Ee_stol,diff_scat_atoms[Ee_stol],count[Ee],"N"};\
}' |\
sort --key=1,2 --key=2n,3 |\
egrep -v "^H" |\
cat >! ${tempfile}diffuse_scatter_atoms.txt
# format: Ee sin(theta)/lambda sum(diff_fraction) N "N"


# now we need the actual atomic form factors
# might as well steal them from CCP4
echo -n >! ${tempfile}Ee_F_compton.txt
foreach atom ( `awk '{print $1}' ${tempfile}diffuse_scatter_atoms.txt | sort -u` )

    # try different variations
    set Ee = "$atom"
    set test = `egrep "^$Ee " ${CLIBD}/atomsf.lib | wc -l`
    if("$test" == "0") then
	set Ee = `echo $atom | awk '{print substr($1,1,1)tolower(substr($1,2,1))}'`
    endif
    set test = `egrep "^$Ee " ${CLIBD}/atomsf.lib | wc -l`
    if("$test" == "0") then
	set Ee = `egrep -i "^$atom" ${CLIBD}/atomsf.lib | head -1 | awk '{print $1}'`
    endif
    set test = `egrep "^$Ee " ${CLIBD}/atomsf.lib | wc -l`
    if("$test" == "0") then
	echo "WARNING: unable to find form factor for $atom "
	continue
    endif

    if("$Ee" != "$atom") echo "using $Ee form-factor for $atom"

    # extract the form factor in 0.01 steps of sin(theta)/lambda
    # also calculate the Compton inline here
    cat ${CLIBD}/atomsf.lib |\
    awk -v Ee="$Ee" 'NF==1 && $1==Ee{getline;Z=$2;c=$3;\
        getline; for(i=1;i<=4;++i) a[i]=$i;\
        getline; for(i=1;i<=4;++i) b[i]=$i;\
        print "Z=",Z;\
        for(s=0;s<=2;s+=0.01){\
        print s,a[1]*exp(-b[1]*s*s)+a[2]*exp(-b[2]*s*s)+a[3]*exp(-b[3]*s*s)+a[4]*exp(-b[4]*s*s)+c;\
    }}' |\
    awk -v atom=$atom -v energy=$energy 'BEGIN{lambda=12398.4245/energy}\
       /^Z=/{Z=$NF;next}\
       {stol=$1;F=$2;\
        sintheta=stol*lambda;if(sintheta>1) next\
        costheta=sqrt(1-sintheta^2);\
        kappa=energy/511000;\
        polar=(1+costheta^2)/2;\
        compton=(Z-F)*polar/(1+kappa*(1-costheta))^2;\
        if(compton<0)compton=0;\
        print atom,stol,F,compton,"COMPTON";\
    }' >> ${tempfile}Ee_F_compton.txt
    # format: Ee sintheta/lambda F compton_electrons

end

# now combine the diffuse scatter fraction sum and atom count with the
# structure factor and Compton "exess electrons"
# diffuse scatter is the difference between the scattering from the atom in the gas phase
# and the scattering from the atom in the crystal lattice
# Compton electrons scatter individually, so this term is NOT squared
cat ${tempfile}diffuse_scatter_atoms.txt ${tempfile}Ee_F_compton.txt |\
awk '{Ee=$1;stol=$2} $NF=="N"{diff[Ee,stol]=$3;N[Ee]=$4;next}\
    $NF=="COMPTON"{F=$3;compton=$4;print Ee,stol,F*F*diff[Ee,stol],compton*N[Ee]}' |\
awk '{dsum[$2]+=$3;csum[$2]+=$4}\
    END{for(stol in dsum){if(dsum[stol]+0>0 && csum[stol]+0>0){\
	print stol,dsum[stol]+csum[stol],dsum[stol],csum[stol]}}}' |\
sort -n >! ${tempfile}incoherent_scatter.stol
# format: stol diffuse(events/sterad/cell) Compton(events/sterad/cell)


# fluence*(re^2)*"molecules"
set I_over_Fsq = `echo $noisy_thomson $ASUs | awk '{print $1*$2}'`

# generate the form factor for the incoherent scattering
set Rmax = `echo $dimx $dimy | awk '{print 2*int(sqrt($1*$1+$2*$2))}'`
set samples = `echo $Rmax | awk '{print $1+1}'`

gnuplot -persist << EOF >&! ${tempfile}incoherent_scatter.log
scale = ${I_over_Fsq}
dist = $distance
pixelsize = $pixel
lambda = $wavelength
mu_air = $mu_air
Rmax   = $Rmax
# shorthand...
c(x)=column(x)
safexp(x) = (x>700 ? 1e300 : (x<-700 ? 0 : exp(x)))

# Bragg angle
theta(stol) = (stol*lambda<=1&&stol>=0?atan2(stol*lambda,sqrt(1-(stol*lambda)**2)):1/0)
# take-off angle from sample
twotheta(stol) = 2*theta(stol)
# path from sample to detector surface (in mm)
airpath(stol) = dist/cos(twotheta(stol))
# size of a pixel in Steradians
pixelSR(stol) = (pixelsize/airpath(stol))**2*cos(2*theta(stol))
# distance of a pixel from the beam center (in pixel units)
pixrad(stol) = airpath(stol)*sin(2*theta(stol))/pixelsize
# intensity in photons/pixel for this geometry
I(stol,fsqr) = scale*safexp(-airpath(stol)/mu_air)*pixelSR(stol)*fsqr
# sanity check
sanepixrad(stol) = (pixrad(stol)<Rmax && twotheta(stol)<pi/2?pixrad(stol):1/0)
saneI(stol,fsqr) = (twotheta(stol)<pi/2?I(stol,fsqr):1/0)
set xlabel "radius (pixels)"
set ylabel "photons/pixel"
set table "${tempfile}plotout"
set output "${tempfile}plotout"
$gnuploterminal
set samples $samples
plot [0:Rmax] '${tempfile}incoherent_scatter.stol' using (sanepixrad(c(1))):(I(c(1),c(2)))
set output
! awk '\$3=="i"' ${tempfile}plotout > ${tempfile}replot
! echo $Rmax | awk '{print 10*\$1,0}' >> ${tempfile}replot
set output "${tempfile}plotout"
plot [0:Rmax] '${tempfile}replot' using 1:2 smooth csplines
EOF

# convert to a format that FIT2D can read in
#peter - edited this to make it work
tail -n+2 ${tempfile}plotout | awk '$3=="i"{print int($1),$2}' |\
awk 'BEGIN{clip=1} $2<=0{clip=0} {print $1,$2*clip}' >! ${tempfile}photons_per_pixel_from_incoherent_scatter.txt

cat << EOF >! ${tempfile}.incoherent_scatter
STORE
INPUT
1-D ASCII
${tempfile}photons_per_pixel_from_incoherent_scatter.txt
VERT
0
0

1
2
EXCHANGE
SYMM
EOF
echo "$beam_center $pixel" |\
 awk '{print $2/$3; print $1/$3}' |\
cat >> ${tempfile}.incoherent_scatter

# don't bother if nothing is there
if(! -s ${tempfile}photons_per_pixel_from_incoherent_scatter.txt) then
    echo "WARNING: incoherent scattering unavailable.  requires gnuplot 4.6 and refined.pdb file in cwd"
    echo -n "" >! ${tempfile}.incoherent_scatter
endif









# non-crystal scattering (water "thickness" should include the water in solvent channels! )
echo -n "" >! ${tempfile}.nonxtal_scatter
set nbg = 0
foreach material ( $nonxtal_materials )
    @ nbg = ( $nbg + 1 )
    if($#nonxtal_thickness < $nbg) then
	echo "ERROR: no thickness specified for $material"
	echo "       unable to add $material scatter"
	continue
    endif
    if($#nonxtal_densities < $nbg) then
	echo "ERROR: no density specified for $material"
	echo "       unable to add $material scatter"
	continue
    endif
    if($#nonxtal_molweight < $nbg) then
	echo "ERROR: no molecular/formula weight specified for $material"
	echo "       unable to add $material scatter"
	continue
    endif
    set thick = `echo "$nonxtal_thickness[$nbg]" | awk '$1+0>1e-3{print $1+0}'`
    set density = `echo "$nonxtal_densities[$nbg]" | awk '$1+0>0{print $1+0}'`
    set molweight = `echo "$nonxtal_molweight[$nbg]" | awk '$1+0>0.1{print $1+0}'`

    if("$thick" == "" || "$density" == "") then
	continue
    endif
    if("$molweight" == "") then
	echo "ERROR: molecular weight of $material cannot be $nonxtal_molweight[$nbg] atomic mass units "
	continue
    endif
    set dir = "."
    if(! -e ${dir}/${material}.stol) then
	set dir = `dirname $0`
    endif
    if(! -e ${dir}/${material}.stol) then
	echo "ERROR: ${material}.stol file does not exist! "
	echo "       unable to add $material scatter"
	continue
    endif

    # calculate the amount of this material in the beam (in MKS units)
    #set volume = `echo $thick $beam_size | awk '{print ($1/1e6)*($2/2/1e6)^2*3.1415}'`
    set volume = `echo $thick $beam_size | awk '{print ($1/1e6)*($2/1e6)^2}'`
    set molecules = `echo 6.022e23 1.66054e-27kg $volume $molweight $density | awk '{NA=$1;amu=$2;V=$3;MW=$4/1000;density=$5*1000;print NA*(density*V/MW)}'`

    # fluence*(re^2)*molecules
    set I_over_Fsq = `echo $noisy_thomson $molecules | awk '{print $1*$2}'`

    # generate the form factor for this material
    # (from a file containing "electrons/molecule" vs sin(theta)/lambda)
    set Rmax = `echo $dimx $dimy | awk '{print 2*int(sqrt($1*$1+$2*$2))}'`
    set samples = `echo $Rmax | awk '{print $1+1}'`
gnuplot -persist << EOF >& ${tempfile}_gnuplot_${material}.log
scale = ${I_over_Fsq}
dist = $distance
pixelsize = $pixel
lambda = $wavelength
mu_air = $mu_air
Rmax = $Rmax
# shorthand...
c(x)=column(x)
safexp(x) = (x>700 ? 1e300 : (x<-700 ? 0 : exp(x)))

# Bragg angle
theta(stol) = (stol*lambda<1&&stol>=0?atan2(stol*lambda,sqrt(1-(stol*lambda)**2)):1/0)
# take-off angle from sample
twotheta(stol) = 2*theta(stol)
# path from sample to detector surface (in mm)
airpath(stol) = dist/cos(twotheta(stol))
# size of a pixel in Steradians
pixelSR(stol) = (pixelsize/airpath(stol))**2*cos(2*theta(stol))
# distance of a pixel from the beam center (in pixel units)
pixrad(stol) = airpath(stol)*sin(2*theta(stol))/pixelsize
# intensity in photons/pixel for this geometry
I(stol,electrons) = scale*safexp(-airpath(stol)/mu_air)*pixelSR(stol)*(electrons)**2
# sanity check
sanepixrad(stol) = (pixrad(stol)<Rmax && twotheta(stol)<pi/2?pixrad(stol):1/0)
saneI(stol,electrons) = (twotheta(stol)<pi/2?I(stol,electrons):1/0)
set xlabel "radius (pixels)"
set ylabel "photons/pixel"
set output "${tempfile}plotout"
set table "${tempfile}plotout"
$gnuploterminal
set samples $samples
plot [0:Rmax] '${dir}/${material}.stol' using (sanepixrad(c(1))):(I(c(1),c(2)))
set output
! awk '\$3=="i"' ${tempfile}plotout > ${tempfile}replot
! echo $Rmax | awk '{print 10*\$1,0}' >> ${tempfile}replot
set output "${tempfile}plotout"
plot [0:Rmax] '${tempfile}replot' using 1:2 smooth csplines
#plot [0:Rmax] '${dir}/${material}.stol' using (pixrad(c(1))):(saneI(c(1),c(2))) smooth csplines
EOF

    # convert to a format that FIT2D can read in
    #peter - edited this to make it work
    tail -n+2 ${tempfile}plotout | awk '$3=="i"{print int($1),$2}' |\
    awk 'BEGIN{clip=1} $2<=0{clip=0} {print $1,$2*clip}' >! ${tempfile}photons_per_pixel_from_${material}.txt


    cat << EOF >> ${tempfile}.nonxtal_scatter
STORE
INPUT
1-D ASCII
${tempfile}photons_per_pixel_from_${material}.txt
VERT
0
0

1
2
EXCHANGE
SYMM
EOF
    echo "$beam_center $pixel" |\
     awk '{print $2/$3; print $1/$3}' |\
    cat >> ${tempfile}.nonxtal_scatter

end









# make ice rings (using relative intensities recorded by Bragg (1922))
echo -n >! ${tempfile}.icerings
foreach icering ( 3.915_1 3.671_10 3.453_2 2.675_1.5 2.260_1 2.065_5 1.925_1 1.528_1.5 1.372_2 1.305_0.25 1.268_0.25 1.165_0.5 )

set d        = `echo $icering | awk -F "[_]" '{print $1}'`
set strength = `echo $icering | awk -F "[_]" '{print $2}'`

echo "$d $icering_width $strength $icering_intensity $wavelength $distance $beam_center $pixel" |\
awk 'BEGIN{PI=2*atan2(1,0)} {\
   d=$1;width=$2;height=$3*$4;lambda=$5;dist=$6;xbeam=$7;ybeam=$8;pix=$9;\
   sintheta=lambda/2/d;if(sintheta>1)next;\
   theta=atan2(sintheta,sqrt(1-sintheta^2));\
   print "RING"; print pix*1000; print pix*1000;\
   print ybeam/pix; print xbeam/pix; print dist;\
   print 0; print 0; print 2*theta/PI*180; print height; print width;\
   print "RING"; print pix*1000; print pix*1000;\
   print ybeam/pix; print xbeam/pix; print dist;\
   print 0; print 0; print 2*theta/PI*180; print height/20; print width*3.5}' |\
cat >> ${tempfile}.icerings

end



if("$icering_intensity" == 0) then
    echo -n "" >! ${tempfile}.icerings
endif
if ( $?NO_NOISE ) then
    echo -n "" >! ${tempfile}.noise
    echo -n "" >! ${tempfile}.ripplenoise
#    echo -n "" >! ${tempfile}.psf
endif
if($?NO_SPOTS) then
    echo -n "" >! ${tempfile}_${batch}.spots
endif
if($?NO_BACKGROUND) then
    echo -n "" >! ${tempfile}.incoherent_scatter
    echo -n "" >! ${tempfile}.nonxtal_scatter
    echo -n "" >! ${tempfile}.saxs
endif

# shame we can't do this without graphics
cat << EOF >! ${tempfile}.polar
polarization
EOF

# generate the full FIT2D script
cat ${tempfile}.start \
${tempfile}_${batch}.spots \
${tempfile}.ripplenoise \
${tempfile}.incoherent_scatter \
${tempfile}.nonxtal_scatter \
${tempfile}.saxs ${tempfile}.icerings \
${tempfile}.beamstop \
${tempfile}.directbeam \
${tempfile}.noise ${tempfile}.psf ${tempfile}.adc \
${tempfile}.windowpane ${tempfile}.end |\
cat >! fit2d_${num}.in

echo "fit2d run..."
cat fit2d_${num}.in |\
fit2d -nographics >! ${tempfile}.fit2d.log

# go easy on NFS
mv ${tempfile}.fit2d.log fit2d_${num}.log

cat ${tempfile}.header |\
awk 'BEGIN{pad=512} {print;pad-=length($0)+1}\
       END{printf "\f"; while(pad>1) {printf " "; --pad}}' |\
cat >! $outfile
if("$status") then
    echo "ERROR: cannot create $outfile"
    goto exit
endif

cat ${tempfile}fit2d$$.bin >> $outfile
rm -f ${tempfile}fit2d$$.bin

msdate.com $start_time
echo "s to generate $outfile"
set start_time = `msdate.com | awk '{print $NF}'`
echo ""

end








########################################
exit:

if($?BAD) then
    echo "ERROR: $BAD"
    exit 9
endif


exit







################################################################################################
Setup:
################################################################################################


set dosemode = 1
set binmode  = 2
set newdark  = 1
set inverse_beam = 0

set name         = ""
set directory    = ""
set runnum       = 1

set firstbatch   = 1
set wedge        = ""
set skip_step    = 1
set skip_start   = 0

set binmode      = 1

set Amatrix = none

set defaultenergy
unset userframes


if(-e mlfsom.params) then
    echo "reading mlfsom.params"
    source mlfsom.params
endif

set energy       = `echo $wavelength | awk '{print 12398.4245/$1}'`
set energies     = ( $energy )

foreach arg ( $* )
    if("$arg" =~ *.mtz) then
	if(-e "$arg") then
	    set mtzfile = $arg
	else
	    echo "WARNING: $arg does not exist."
	endif
    endif
    if("$arg" =~ *.pdb) then
	if(-e "$arg") then
	    set pdbfile = $arg
	else
	    echo "WARNING: $arg does not exist."
	endif
    endif
    if("$arg" =~ *.img) then
	set directory = `dirname $arg`
	if (-e "$directory") set directory = `cd $directory ; pwd`

	set name = `basename $arg .img`
	set num  = `echo $name | awk -F "_" '{print $NF}'`
	set firstbatch    = `echo $num  | awk '{print $NF+0}'`
	set name = `basename $name _$num`
#	set name = `basename $name _0`
#	set name = `basename $name _1`
#	set name = `echo $name | awk -F "_" '{for(i=1;i<=NF-2;++i) {printf "%s", $i; if(i<NF-2) printf "_"}}'`
    endif
    if("$arg" == "nodark") set newdark = 0
    if("$arg" == "newdark") set newdark = 1
    if("$arg" == "dark") set newdark = 1
    if("$arg" == "cache") set CACHE
    if("$arg" == "nocache") unset CACHE

    if("$arg" == "dose") set dosemode = 1
    if("$arg" == "dosemode") set dosemode = 1
    if("$arg" == "-dosemode") set dosemode = 1
    if("$arg" == "nodose") set dosemode = 0
    if("$arg" == "nodosemode") set dosemode = 0
    if("$arg" == "-nodosemode") set dosemode = 0

    set number = `echo $arg | awk '$1+0>0.001{print $1+0}'`
    if("$arg" =~ *deg && "$number" != "") set arg = "phi=$number"
    if("$arg" =~ *mm && "$number" != "") set arg = "dist=$number"
    if("$arg" =~ *s && "$number" != "")  set arg = "time=$number"
    if("$arg" =~ *A && "$number" != "")  set arg = "wave=$number"
    if("$arg" =~ *eV && "$number" != "") set arg = "energy=$number"

    set afterequals = `echo $arg | awk -F "=" '{print $NF}'`
    set afternum = `echo $arg | awk -F "=" '$NF+0>0.0001{print $NF+0}'`
    set afternums = `echo $arg | awk -F "=" '$NF ~ /^[0-9-]/{print $NF}' | awk -F "," '{print $1,$2,$3}'`

    if("$arg" =~ firstbatch* && "$afternum" != "") then
	set temp = `echo $afternum | awk '$1>=0{print int($1)}'`
	if("$temp" != "") set firstbatch = $afternum
    endif

    if("$arg" =~ mosa* && "$afternum" != "") then
	set temp = `echo $afternum | awk '$1>0 && $1<10'`
	if("$temp" != "") set mosaic = $afternum
    endif
    if("$arg" =~ rock* && "$afterequals" != "") then
	set temp = `echo $afternum | awk '$1>0 && $1<10'`
	if("$afterequals" =~ tan*) set rocking_curve = tanh
	if("$afterequals" =~ gau*) set rocking_curve = gauss
	if("$afterequals" =~ norm*) set rocking_curve = gauss
	if("$afterequals" =~ lore*) set rocking_curve = lorentz
	if("$afterequals" =~ squa*) set rocking_curve = square
	if("$afterequals" =~ toph*) set rocking_curve = tophat
	if("$afterequals" =~ dis[ck]) set rocking_curve = disk
    endif
    if("$arg" =~ halflife* && "$afternum" != "") then
	set temp = `echo $afternum | awk '$1>0'`
	if("$temp" != "") set decay_timeconstant = $afternum
    endif
    if("$arg" =~ spotw* && "$afternum" != "") then
	set temp = `echo $afternum | awk '$1>0 && $1<10'`
	if("$temp" != "") set spotwidth = $afternum
    endif
    if("$arg" =~ xtal[_s]* && "$afternum" != "") then
	set temp = `echo $afternum | awk '$1>=0'`
	if("$temp" != "") set xtal_size = $afternum
    endif
    if("$arg" =~ beam[_s]* && "$afternum" != "") then
	set temp = `echo $afternum | awk '$1>=0'`
	if("$temp" != "") set beam_size = $afternum
    endif
    if("$arg" =~ drop[_s]* && "$afternum" != "") then
	set temp = `echo $afternum | awk '$1>=0'`
	if("$temp" != "") set drop_size = $afternum
    endif
    if("$arg" =~ icering_w* && "$afternum" != "") then
	set temp = `echo $afternum | awk '$1>=0'`
	if("$temp" != "") set icering_width = $afternum
	continue
    endif
    if("$arg" =~ ice* && "$afternum" != "") then
	set temp = `echo $afternum | awk '$1>=0'`
	if("$temp" != "") set icering_intensity = $afternum
	continue
    endif
    if("$arg" =~ misset* && "$arg" =~ *[0-9]*) then
	if ($#afternums == 3) set missets = ( $afternums )
    endif
    if("$arg" =~ matrix* && "$afterequals" =~ *.mat || "$arg" =~ *.mat) then
	set matfile = "$afterequals"
	if(! -e "$matfile" && "$arg" =~ *.mat) then
            set matfile = "$arg"
	endif
	if(-e "$matfile") then
	    set Amatrix = `head -3 $matfile`
	    set rcell   = `head -3 $matfile | awk '{++j;for(i=1;i<=3;++i){A[i,j]=$i}} END{for(i=1;i<=3;++i)printf "%.8f ", sqrt(A[i,1]^2+A[i,2]^2+A[i,3]^2); print ""}'`
	    awk 'NF==6' $matfile | head -1 |\
	    awk '{PI=2*atan2(1,0);\
		    a=$1;b=$2;c=$3;alpha=$4*PI/180;beta=$5*PI/180;gamma=$6*PI/180;\
		    s=(alpha+beta+gamma)/2;skew=sin(s)*sin(s-alpha)*sin(s-beta)*sin(s-gamma);\
		    if(skew<0) skew=-skew; if(skew==0) skew=0.001;\
		    Volume = 2*a*b*c*sqrt(skew);\
		    a_star = b*c*sin(alpha)/Volume;\
		    b_star = c*a*sin(beta)/Volume;\
		    c_star = a*b*sin(gamma)/Volume;\
		    print a_star, b_star, c_star}' |\
	    cat >! ${tempfile}rcell
	    set true_rcell = `cat ${tempfile}rcell`
	    rm -f ${tempfile}rcell
	    set old_wavelength = `echo $rcell $true_rcell | awk '{print ($1/$4+$2/$5+$3/$6)/3}'`
	    set Amatrix = `head -3 $matfile | awk -v l=$old_wavelength '{print $1/l,$2/l,$3/l}'`

	    set missets = ""
	else
	    echo "WARNING: $matfile does not exist"
	endif
    endif
    if("$arg" =~ beam* && "$arg" =~ *[0-9]*) then
	if ($#afternums == 2) set beam_center = ( $afternums )
    endif

    if("$arg" =~ id=* && "$arg" =~ *[0-9]*) then
      set id = `echo $arg | awk -F "=" '{print $NF+0}'`
      set tempfile = ${CCP4_SCR}/mlfsom_tempfile${id}
    endif
    if("$arg" =~ phi=* && "$arg" =~ *[0-9]*) then
	set phi0 = `echo $arg | awk -F "=" '{print $NF+0}'`
    endif
    if("$arg" =~ start=* && "$arg" =~ *[0-9]*) then
	set phi0 = `echo $arg | awk -F "=" '{print $NF+0}'`
    endif
    if("$arg" =~ time=* && "$afternum" != "") then
	set temp = `echo $afternum | awk '$1>0'`
	set exposure = $afternum
    endif
   if("$arg" =~ expos* && "$afternum" != "") then
	set temp = `echo $afternum | awk '$1>0'`
	set exposure = $afternum
    endif
    if(("$arg" =~ delta=* || "$arg" =~ osc=*) && "$afternum" != "") then
	set temp = `echo $afternum | awk '$1>0 && $1<361'`
	if("$temp" != "") set osc = $afternum
    endif
    if(("$arg" =~ spindle_swing=* || "$arg" =~ spindle_swing=*) && "$afternum" != "") then
        set temp = `echo $afternum | awk '$1>0 && $1<361'`
        if("$temp" != "") set spindle_swing = $afternum
    endif
    if("$arg" =~ wedge=* && "$afternum" != "") then
	set temp = `echo $afternum | awk '$1>0 && $1<361'`
	if("$temp" != "") set wedge = $afternum
    endif
    if("$arg" =~ range=* && "$afternum" != "") then
	set temp = `echo $afternum | awk '$1>0.01 && $1<361'`
	if("$temp" != "") set range = $afternum
	unset userframes
    endif
    if("$arg" =~ frames=* && "$afternum" != "") then
	set temp = `echo $afternum | awk '$1>0.01 && $1<361'`
	if("$temp" != "") set userframes = $afternum
    endif
    if("$arg" =~ skip_step=* && "$afternum" != "") then
	set temp = `echo $afternum | awk '$1==int($1) && $1>=1'`
	if("$temp" != "") set skip_step = $afternum
    endif
    if("$arg" =~ skip_start=* && "$afternum" != "") then
	set temp = `echo $afternum | awk '$1==int($1) && $1>=1'`
	if("$temp" != "") set skip_start = $afternum
    endif


    if("$arg" =~ dist* && "$afternum" != "") then
	set temp = `echo $afternum | awk '$1>1 && $1<20000'`
	if("$temp" != "") set distance = $afternum
    endif
    if("$arg" =~ energ* && "$afternum" != "") then
	set temp = `echo $afternum | awk '$1>1000 && $1<2000000'`
	if("$temp" != "") then
	    if($?defaultenergy) set energies = ""
	    unset defaultenergy
	    set energies = ( $energies $afternum )
	endif
    endif
    if("$arg" =~ wave* && "$afternum" != "") then
	set temp = `echo $afternum | awk '$1>0.001 && $1<50{printf "%.2f", 12398.4245/$1}'`
	if("$temp" != "") then
	    if($?defaultenergy) set energies = ""
	    unset defaultenergy
	    set energies = ( $energies $temp )
	else
	    echo "WARNING: ignored unreasonable wavelength: $afternum"
	endif
    endif
    if("$arg" =~ reso* && "$afternum" != "") then
	set temp = `echo $afternum | awk '$1>0.1 && $1<500{printf "%f", $1}'`
	if("$temp" != "") then
	    set USER_res = "$temp"
	endif
    endif
    if("$arg" =~ flux* && "$afternum" != "") then
	set temp = `echo $afternum | awk '$1>0 {printf "%.3g", $1}'`
	if("$temp" != "") then
	    set flux = "$temp"
	endif
    endif
    if("$arg" =~ -inverse*) then
	set inverse_beam = 1
    endif
    if("$arg" =~ -noinverse*) then
	set inverse_beam = 0
    endif
end

# count wavelengths
set energies = ( $energies )
set nwaves   = $#energies

# a few sanity checks
if("$name" == "" && "$outprefix" != "") set name = `basename $outprefix`
if("$directory" == "" && "$outprefix" != "") set directory = `dirname $outprefix`
if("$name" == "") set name = fake
if("$directory" == "") set directory = /data/${user}
if("$firstbatch" == "") set firstbatch = 1

# check that directory is realistic
if(! -w $directory && ! $?NORUN) then
    mkdir -p $directory >& /dev/null
endif
if(! -w $directory && ! $?NORUN) then
    set test = `echo $directory | awk -F "/" '{print "/"$2"/"$3}'`
    if(! -w "$test") then
	# this is not going to work, replace second directory with username
	set directory = `echo $directory | awk -v user=$user -F "/" '{$3=user ; for(d=2;d<=NF;++d) printf "/%s",$d}'`
    endif
endif

set num = `echo $firstbatch | awk '{printf "%03d", $1+0}'`

if($?userframes) then
    set range = `echo $userframes $osc | awk '{print $1*$2}'`
    unset userframes
endif

# calculate the ending phi value
set phiend = `echo $phi0 $range | awk '{print $1+$2}'`
set phiend  = `echo $phi0 $phiend $osc | awk '{n=($2-$1)/$3} n+0!=int(n)+0{n=int(n+1)} {print $1+n*$3}'`
set frames  = `echo $phi0 $phiend $osc 0 $nwaves $inverse_beam | awk '{n=($2-$1)/$3 - $4; print n*$5*($6+1)}'`
set firstnum = `echo $firstbatch | awk '{printf "%03d",$1}'`
set lastnum = `echo $firstbatch $frames | awk '{printf "%03d",$1+$2-1}'`

if( `echo $osc $wedge | awk '{print ($1>$2)}'` ) set wedge = $osc

set firstimage = ${directory}/${name}_$runnum
set lastimage  = ${directory}/${name}_$runnum
if($#energies > 1) then
    set firstimage = "${firstimage}_E1"
    set lastimage  = "${lastimage}_E$#energies"
endif
set firstimage = "${firstimage}_${firstnum}.img"
set lastimage = "${lastimage}_${lastnum}.img"

set energies = ( $energies )
set nwaves   = $#energies
# pad up the energies
while ( $#energies < 5 )
    set energies = ( $energies 0 )
end
set energy = $energies[1]
if( "$energy" != "" ) set wavelength = `echo $energy | awk '{print 12398.4245/$1}'`

if(! $?brightness) set brightness
if("$brightness" == "") then
    set brightness = `echo $flux $beam_size | awk '{print $1/$2/$2}'`
endif

# calculate (electrons/photon) / (electrons/ADU) =  ADU/photon
set gain = `echo $EOgain $ampgain | awk '{print $1/$2}'`


# remove units so they don't mess up fit2d
set dimx = `echo $dimx | awk '{print $1+0}'`
set dimy = `echo $dimy | awk '{print $1+0}'`
set pixel = `echo $pixel | awk '{print $1+0}'`
set ampgain = `echo $ampgain | awk '{print $1+0}'`
set EOgain = `echo $EOgain | awk '{print $1+0}'`
set gain = `echo $gain | awk '{print $1+0}'`
set readnoise = `echo $readnoise | awk '{print $1+0}'`
set adcoffset = `echo $adcoffset | awk '{print $1+0}'`
set psf = `echo $psf | awk '{print $1+0}'`
set readnoise = `echo $readnoise | awk '{print $1+0}'`
set spotwidth = `echo $spotwidth | awk '{print $1+0}'`
set exposure = `echo $exposure | awk '{print $1+0}'`
set brightness = `echo $brightness | awk '{print $1+0}'`
set xtal_size = `echo $xtal_size | awk '{print $1+0}'`
set beam_size = `echo $beam_size | awk '{print $1+0}'`
set wavelength = `echo $wavelength | awk '{print $1+0}'`
set dispersion = `echo $dispersion | awk '{print $1+0}'`
set polarization = `echo $polarization | awk '{print $1+0}'`
set distance = `echo $distance | awk '{print $1+0}'`
set twotheta = `echo $twotheta | awk '{print $1+0}'`
set mosaic = `echo $mosaic | awk '{print $1+0}'`
set phi0 = `echo $phi0 | awk '{print $1+0}'`
set osc = `echo $osc | awk '{print $1+0}'`
set spindle_swing = `echo $spindle_swing | awk '{print $1+0}'`
set frames = `echo $frames | awk '{print $1+0}'`
#set airscatter_path = `echo $airscatter_path | awk '{print $1+0}'`
set icering_intensity = `echo $icering_intensity | awk '{print $1+0}'`
set SAXS_intensity = `echo $SAXS_intensity | awk '{print $1+0}'`

set mu_air = `echo $mu_air | awk '{print $1+0}'`

# do a few conversions
set adxv_center = `echo $beam_center $pixel $dimx $dimy | awk '{w=$5*$3; print $2+0,w-$1}'`

set outprefix = "$directory/$name"

goto Return_from_Setup
