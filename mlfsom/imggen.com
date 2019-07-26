#! /bin/tcsh -f
#
#	render images by running mlfsom.com on multiple CPUs
#
#	requires mlfsom.params and ccp4.setup (for default shell) in the current working directory

cat << EOF >! ccp4.setup
source /programs/ccp4-6.1.3/ccp4-6.1.3/include/ccp4.setup-csh
EOF

set tempfile
awk '/^set tempfile/{gsub("\\$\\$","");print}' mlfsom.com >! mlfsom.sourceme
source mlfsom.sourceme
rm -f mlfsom.sourceme
if("$tempfile" == "") then
    echo "error: mlfsom.com corrupt"
    exit 9
endif

set pwd = `pwd`


alias rsh ssh -x

source mlfsom.params
set outdir = `dirname $outprefix`
set frames = `echo $range $osc | awk '{print int($1/$2)}'`

# kill old jobs
echo "killing any old jobs..."
set machines = `cat machines.txt`
foreach machine ( $machines )
    rsh -n $machine "killall -g mlfsom.com ; rm -f ${tempfile}*"
end
echo "old jobs killed"
sleep 3
echo "starting new ones..."

rm -rf ${outdir}
set fb = 1
find ${outdir} -size -5 -exec rm -f \{\} \;

rm -f ${tempfile}*rip*
if(! -e ./ripplemask.f2d) then
    rm -f ./ripplemask.f2d
    mkdir -p ${outdir}
    touch ${outprefix}_001.img
    setenv CACHE images
    ./mlfsom.com frames=1 firstbatch=1
    echo "saving the ripple mask as ./ripplemask.f2d"
    mv ${tempfile}*_ripplemask.f2d ./ripplemask.f2d
    rm -f ${outprefix}_001.img
endif

set machines = `cat machines.txt`
set i = 0
foreach machine ( $machines )

@ i = ( $i + 1 )

if($fb > $frames) continue

set start = `echo $phi0 $fb $osc | awk '{print $1+($2-1)*$3}'`
set toend = `echo $frames $fb | awk '{print $1-$2+1}'`
if( $toend < 1 ) set toend = 1

#rsh -n $machine "rm -f ${tempfile}*"
echo "$machine  firstbatch=$fb start=$start"
rsh -n $machine "cd $pwd ; source ccp4.setup ; ./mlfsom.com frames=$toend firstbatch=$fb start=$start" >&! ${machine}_${i}.log &
@ fb = ( $fb + 1 )
sleep 1

end

wait

# fill in anything that got messed up?
set filesize = `echo $dim | awk '{print $1*$1*2+512}'`
find ${outdir} -size -${filesize}c -type f -exec rm -f \{\} \;
./mlfsom.com frames=$frames







