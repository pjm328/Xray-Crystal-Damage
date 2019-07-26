#! /usr/bin/awk -f
#
#	hkl2XYphi.awk				-James Holton  1-30-12
#
#	take a crystal orientation and camera parameters
#	and predict where the spots are going to fall
#
#	
#	z is the rotation axis
#
# example input:
# CELL 74 74 34 90 90 90
# MISSET 10.0 20.0 30.0
# RESOL 1.8
# DIST 120
# WAVE 1.00
# PIXEL 0.1024
# POLAR 0.9
# BEAM  152.5 152.5
#


BEGIN{
PI = 2*atan2(1,0)
RTD = 180/PI

AX = 0.0147471  
AY = 0.000000  
AZ = 0.000000  
BX = 0.00737357
BY = 0.0127714
BZ = 0.000000
CX = 0.000000
CY = 0.000000
CZ = 0.00411184

# default rotation matrix
uxx=1;uxy=0;uxz=0; 
uyx=0;uyy=1;uyz=0;
uzx=0;uzy=0;uzz=1;

d_min = 2
dist = 100
lambda0 = 1
width = 2048*0.1024
xbeam = ybeam = width/2
pixel = 0.1024

# typical values for 12660 eV
mu_air       = 3220
# expected values for ADSC: Gd2O2S:Tb phospor and Al/plastic window
window_thick = 12.7
mu_window    = 610
det_thick    = 40
mu_det       = 10.9
muen_det     = 11.1
# mean-free path of a visible photon in the phosphor
phosphor_mfp = 100

# minimum number of pixels in a "feature"
speedup_pix = 1.5
# number of kernels within mosaic spread
speedup_phi = 1

max_subspots = 1000

# spindle non-orthogonality
spindle_swing = 0.0

# anisotropic mosaic spread?
mosaic_tanmult = 1.0
phimosaic_divisions = 0
tanmosaic_divisions = 0
}

/^CPU/{thiscpu=$2+0;CPUs=$3+0}

/^RESOL/{d_min=$2}
/^DIST/{dist=$2}
/^TWOTHETA/{detector_twotheta=-$2/RTD}
/^TILT/{detector_tilt=$2/RTD}
/^TWIST/{detector_twist=-$2/RTD}
/^CCOMEGA/{detector_omega=$2/RTD}
/^WAVE/{lambda0=$2}
/^ENERGY/{lambda0=12398.42/$2}
/^DISPER/{dispersion=$2;correlated_dispersion=$3}
/^DIVE/{Hdiv=$2/RTD;Vdiv=$3/RTD}
#/^MOSA/{mosaic_spread=$2/RTD;if($3+0>0 && $2+0>0){mosaic_tanmult=$3/$2}}
/^MOSA/{mosaic_spread=$2/RTD}
/^MOSDIV/{phimosaic_divisions=$2;tanmosaic_divisions=$3+0}
/^DISPDIV/{disp_divisions=$2}
/^DIVDIV/{hdiv_divisions=$2;vdiv_divisions=$3}
/^PIX/{pixel=$2}
/^POLAR/{polar_fac=$NF}
/^BEAM /{xbeam=$2;ybeam=$3}
/^BEAMSIZE/{beam_size=$2/1000}
/^XTALSIZE/{xtal_size=$2/1000}
/^DETSIZE/{det_size=$2}
/^HKL/{user_hkl=1;user_h=$2;user_k=$3;user_l=$4}

/^MAX_SUB_SPOTS/{max_subspots=$2}
/^SPEEDUP /{speedup_phi=$2;speedup_pix=$2}
/^SPEEDUP_PIX/{speedup_pix=$2}
/^SPEEDUP_PHI/{speedup_phi=$2}

/^SPINDLE/{spindle_swing=$2/RTD}

# absorption model
# air attenuation length (mu) should be in mm (because um would just be silly)
/^AIR_MU/{mu_air=$2}
# all other mu and thicknesses in um
/^WINDOW_THICK/{window_thick=$2}
/^WINDOW_MU/{mu_window=$2}
/^DET_THICK/{det_thick=$2}
/^DET_MU /{mu_det=$2}
/^DET_MUEN/{muen_det=$2}
/^DET_MFP/{phosphor_mfp=$2}

/^ABSORB/{++absorbers;
    # depth vector at phi=0 and attenuation coefficient
    absorber_x[absorbers]=$2;
    absorber_y[absorbers]=$3;
    absorber_z[absorbers]=$4;
    if($5+0>0) {
	absorber_mu[absorbers]=$5;
    }else{
	# re-use the last one
	absorber_mu[absorbers]=absorber_mu[absorbers-1];
    }
    # units of mu should be the same as the xyz units
}


# take pseudo-mosflm input
/^MATRIX/ && NF==10{AX=$2;BX=$3;CX=$4;AY=$5;BY=$6;CY=$7;AZ=$8;BZ=$9;CZ=$10; Amatrix=1}

# make a B-matrix out of any cell
/^CELL/ && NF==7{
    a=$2;b=$3;c=$4;alpha=$5/RTD;beta=$6/RTD;gamma=$7/RTD;

    # if cell is defined, create a B-matrix
    s=(alpha+beta+gamma)/2;skew=sin(s)*sin(s-alpha)*sin(s-beta)*sin(s-gamma);
    if(skew<0) skew=-skew; if(skew==0) skew=0.001;
    Volume = 2*a*b*c*sqrt(skew);
    a_star = b*c*sin(alpha)/Volume;
    b_star = c*a*sin(beta)/Volume;
    c_star = a*b*sin(gamma)/Volume;
    cos_alphastar = (cos(beta)*cos(gamma)-cos(alpha))/(sin(beta)*sin(gamma));
    cos_betastar  = (cos(gamma)*cos(alpha)-cos(beta))/(sin(gamma)*sin(alpha));
    cos_gammastar = (cos(alpha)*cos(beta)-cos(gamma))/(sin(alpha)*sin(beta));
    sin_alphastar = sqrt(1-cos_alphastar^2);
    sin_betastar  = sqrt(1-cos_betastar^2);
    sin_gammastar = sqrt(1-cos_gammastar^2);
    B_AX=a_star;
    B_BX=b_star*cos_gammastar;
    B_CX=c_star*cos_betastar;
    B_AY=0;
    B_BY=b_star*sin_gammastar;
    B_CY=c_star*(cos_alphastar-cos_betastar*cos_gammastar)/sin_gammastar;
    B_AZ=0;
    B_BZ=0;
    B_CZ=c_star*Volume/(a*b*c*sin_gammastar);
    Bmatrix=1;
}

# multiply any U-matrix by the current U-matrix
/^UMAT/ && NF==10{
    rxx=$2;rxy=$3;rxz=$4; ryx=$5;ryy=$6;ryz=$7; rzx=$8;rzy=$9;rzz=$10;

    # multiply U = U*R
    oxx=uxx; oxy=uxy; oxz=uxz;
    oyx=uyx; oyy=uyy; oyz=uyz;
    ozx=uzx; ozy=uzy; ozz=uzz;
    uxx = oxx*rxx + oxy*ryx + oxz*rzx;
    uxy = oxx*rxy + oxy*ryy + oxz*rzy;
    uxz = oxx*rxz + oxy*ryz + oxz*rzz;
    uyx = oyx*rxx + oyy*ryx + oyz*rzx;
    uyy = oyx*rxy + oyy*ryy + oyz*rzy;
    uyz = oyx*rxz + oyy*ryz + oyz*rzz;
    uzx = ozx*rxx + ozy*ryx + ozz*rzx;
    uzy = ozx*rxy + ozy*ryy + ozz*rzy;
    uzz = ozx*rxz + ozy*ryz + ozz*rzz;
   
    Umatrix=1;
}

# use misseting angles (in the same way mosflm does)
/^MISSET/ && NF>=4{
    rotX=-$2/RTD; rotY=-$3/RTD; rotZ=-$4/RTD;
    
    # rotate each basis vector in the U matrix
    rotate(uxx,uxy,uxz,rotX,rotY,rotZ);
    uxx=new_x;uxy=new_y;uxz=new_z;
    rotate(uyx,uyy,uyz,rotX,rotY,rotZ);
    uyx=new_x;uyy=new_y;uyz=new_z;
    rotate(uzx,uzy,uzz,rotX,rotY,rotZ);
    uzx=new_x;uzy=new_y;uzz=new_z;

    Umatrix=1;
}


END {

# use B matrix as the A matrix if we have it
if ( Bmatrix && ! Amatrix ) {
AX=B_AX; BX=B_BX; CX=B_CX;
AY=B_AY; BY=B_BY; CY=B_CY;
AZ=B_AZ; BZ=B_BZ; CZ=B_CZ;
Amatrix = 1
}
if(! Amatrix) exit

# rotate A matrix about by the U matrix (or misseting angles)
nxx=AX; nxy=BX; nxz=CX;
nyx=AY; nyy=BY; nyz=CY;
nzx=AZ; nzy=BZ; nzz=CZ;
AX = uxx*nxx + uxy*nyx + uxz*nzx;
BX = uxx*nxy + uxy*nyy + uxz*nzy;
CX = uxx*nxz + uxy*nyz + uxz*nzz;
AY = uyx*nxx + uyy*nyx + uyz*nzx;
BY = uyx*nxy + uyy*nyy + uyz*nzy;
CY = uyx*nxz + uyy*nyz + uyz*nzz;
AZ = uzx*nxx + uzy*nyx + uzz*nzx;
BZ = uzx*nxy + uzy*nyy + uzz*nzy;
CZ = uzx*nxz + uzy*nyz + uzz*nzz;

#    printf "%6.3f %6.3f %6.3f\n", uxx,uxy,uxz
#    printf "%6.3f %6.3f %6.3f\n", uyx,uyy,uyz
#    printf "%6.3f %6.3f %6.3f\n", uzx,uzy,uzz
#    printf "%6.3f %6.3f %6.3f\n", AX+0.00001,BX+0.00001,CX+0.00001
#    printf "%6.3f %6.3f %6.3f\n", AY+0.00001,BY+0.00001,CY+0.00001
#    printf "%6.3f %6.3f %6.3f\n", AZ+0.00001,BZ+0.00001,CZ+0.00001


lambda_star0=1/lambda0
lambda_star0_sqr=lambda_star0^2



detector_tilt = detector_tilt + detector_twotheta;
detector_twist = detector_twist - spindle_swing;

# unit vectors in the direction of each scanning axis of the detector
detector_X_x = 0
detector_X_y = -1
detector_X_z = 0

detector_Y_x = 0
detector_Y_y = 0
detector_Y_z = 1

# normal to the detector face
detector_Z_x = 1
detector_Z_y = 0
detector_Z_z = 0


# rotate the detector around lab frame axes
rotate(detector_X_x,detector_X_y,detector_X_z,-detector_omega,0,0);
detector_X_x=new_x;detector_X_y=new_y;detector_X_z=new_z;
rotate(detector_Y_x,detector_Y_y,detector_Y_z,-detector_omega,0,0);
detector_Y_x=new_x;detector_Y_y=new_y;detector_Y_z=new_z;
rotate(detector_Z_x,detector_Z_y,detector_Z_z,-detector_omega,0,0);
detector_Z_x=new_x;detector_Z_y=new_y;detector_Z_z=new_z;

rotate(detector_X_x,detector_X_y,detector_X_z,0,-detector_twist,-detector_tilt);
detector_X_x=new_x;detector_X_y=new_y;detector_X_z=new_z;
rotate(detector_Y_x,detector_Y_y,detector_Y_z,0,-detector_twist,-detector_tilt);
detector_Y_x=new_x;detector_Y_y=new_y;detector_Y_z=new_z;
rotate(detector_Z_x,detector_Z_y,detector_Z_z,0,-detector_twist,-detector_tilt);
detector_Z_x=new_x;detector_Z_y=new_y;detector_Z_z=new_z;



h_max = int(a/d_min); h_min = -h_max;
k_max = int(b/d_min); k_min = -k_max;
l_max = int(c/d_min); l_min = -l_max;

if(user_hkl) { h_max=h_min=user_h ; k_max=k_min=user_k ; l_max=l_min=user_l; }

# loop over all potentially interesting indicies
for(h=h_min;h<=h_max;h+=1){
for(k=k_min;k<=k_max;k+=1){
for(l=l_min;l<=l_max;l+=1){



# find the point in reciprocal space where this HKL resides at phi=0
x=h*AX + k*BX + l*CX;
y=h*AY + k*BY + l*CY;
z=h*AZ + k*BZ + l*CZ;

# reciprocal-space magnitude of this HKL vector
d_star = sqrt(x^2 + y^2 + z^2)

# skip stuff outside the resolution of interest
if(d_star > 0){
    if(1/d_star < d_min && l>0) break
    if(1/d_star < d_min) continue
}
# skip stuff outside the limiting sphere
if(d_star > 2*lambda_star0) {
    # definitely outside Ewald sphere
    continue
}

# skip every so often if we are running on one of several CPUs
if(CPUs>1) {
    cpu=cpu%CPUs+1;
    if(cpu != thiscpu) continue;
}

# the same spot will appear on either side of rotation axis
# so there will be two intercepts
for(sign=-1;sign<=1;sign+=2) {

if(user_hkl) sign = 1;

# compute the position of the center of the relp
XYphi(h,k,l,0,0,0,0,0);

# check to see if we are in the blind region

if(Xdet=="nan" && required_angle != "nan") {
    # the center is in the blind region
    test_mosaic = required_angle*1.000001;
    if(abs(test_mosaic) < mosaic_spread+Hdiv) {
	XYphi(h,k,l,0,test_mosaic,0,0,0);
    }
}
if(Xdet=="nan" && required_lambda_dev != "nan") {
    # this corner is in the blind region
    test_disp = required_lambda_dev*1.00001;
    if(test_disp > -dispersion && test_disp < dispersion) {
	XYphi(h,k,l,0,0,test_disp,0,0);
    }
}

if(Xdet=="nan") {
    # too far away
#    print "GOTHERE miss",h,k,l,sign,required_angle,required_lambda
    continue
}

if(det_size != "") {
    # don't bother if the ray will miss the detector completely
    if(Xdet < -10) continue
    if(Ydet < -10) continue
    if(Xdet*pixel > det_size+10) continue
    if(Ydet*pixel > det_size+10) continue
}

# too many divisions by zero with this...
if(Lorentz==0) {
    continue
}

# cannot hit the detector from behind!
if(cos(oblique_angle) <= 0) continue



# try to decide how many steps to break each parameter into

# evaluate the X-Y and phi position of this h,k,l,sign at every combination
# of the extreme values of each "width" parameter

spot_corner=0;
max_dpix=max_dphi=max_pix_split=max_phi_split=0;
max_dpix_dphi=0;
max_dpix_tanm=max_dpix_disp=max_dpix_hdiv=max_dpix_vdiv=0;
max_dphi_tanm=max_dphi_disp=max_dphi_hdiv=max_dphi_vdiv=0;
for(stanm=-1;stanm<=1;stanm+=2){
for(sdisp=-1;sdisp<=1;sdisp+=2){
for(shdiv=-1;shdiv<=1;shdiv+=2){
for(svdiv=-1;svdiv<=1;svdiv+=2){

    # refresh parameter values (may get modified below)
    dtanm=stanm*(mosaic_spread+1e-30)/2;
    ddisp=sdisp*(dispersion+1e-30)/2;
    dhdiv=shdiv*(Hdiv+1e-30)/2;
    dvdiv=svdiv*(Vdiv+1e-30)/2;

    # compute the position of this "corner" of the relp
    XYphi(h,k,l,0,dtanm,ddisp,dhdiv,dvdiv);

    # check to see if we are in the blind region

    if(Xdet=="nan" && required_angle != "nan") {
        # this corner is in the blind region
	new_dtanm = dtanm+required_angle*1.000001;
	if(new_dtanm > -mosaic_spread/2 && new_dtanm < mosaic_spread/2) {
	    dtanm = new_dtanm;
	    XYphi(h,k,l,0,dtanm,ddisp,dhdiv,dvdiv);
#	    if(Xdet=="nan") print "GOTHERE: fuck mosaic"
	}
    }

    if(Xdet=="nan" && required_angle != "nan") {
        # this corner is in the blind region
	new_dhdiv = dhdiv-required_angle*1.000001;
	if(new_dhdiv > -Hdiv/2 && new_dhdiv < Hdiv/2) {
	    dhdiv = new_dhdiv;
	    XYphi(h,k,l,0,dtanm,ddisp,dhdiv,dvdiv);
#	    if(Xdet=="nan") print "GOTHERE: fuck hdiv"
	}
    }
    if(Xdet=="nan" && required_lambda_dev != "nan") {
        # this corner is in the blind region
	new_ddisp = required_lambda_dev*1.00001;
	if(new_ddisp > -dispersion/2 && new_ddisp < dispersion/2) {
	    ddisp = new_ddisp;
	    XYphi(h,k,l,0,dtanm,ddisp,dhdiv,dvdiv);
#	    if(Xdet=="nan") print "GOTHERE: fuck disp"
	}
    }

    if(Xdet=="nan"){
	# uhh... now what?
    }

    X_corner[spot_corner]=Xdet;
    Y_corner[spot_corner]=Ydet;
    phi_corner[spot_corner]=phi;
    dtanm_corner[spot_corner]=dtanm;
    ddisp_corner[spot_corner]=ddisp;
    dhdiv_corner[spot_corner]=dhdiv;
    dvdiv_corner[spot_corner]=dvdiv;

    if(Xdet != "nan") {
        for(sc=0;sc<spot_corner;++sc) {
	    if(X_corner[sc]!="nan") {
	        dpix = sqrt((X_corner[sc]-X_corner[spot_corner])^2\
                           +(Y_corner[sc]-Y_corner[spot_corner])^2);
		dphi = sqrt((phi_corner[sc]-phi_corner[spot_corner])^2);
		while(dphi>PI){dphi-=2*PI};
		if(dphi!=0)dpix_dphi=dpix/dphi;
		ddtanm=dtanm_corner[spot_corner]-dtanm_corner[sc];
		dddisp=ddisp_corner[spot_corner]-ddisp_corner[sc];
		ddhdiv=dhdiv_corner[spot_corner]-dhdiv_corner[sc];
		ddvdiv=dvdiv_corner[spot_corner]-dvdiv_corner[sc];
#print "GOTHERE",dpix,dtanm,ddisp,dhdiv,dvdiv;
#print "GOTHERE",dpix,dtanm_corner[sc],ddisp_corner[sc],dhdiv_corner[sc],dvdiv_corner[sc];
#print "GOTHERE1",dpix,dphi,ddtanm,dddisp,ddhdiv,ddvdiv;
	        if(dpix>max_dpix) {
		    max_dpix=dpix;
		    max_sc1=spot_corner;
		    max_sc2=sc;
		}
	        if(dphi>max_dphi) {
		    max_dphi=dphi;
		}
	        if(dpix_dphi>max_dpix_dphi) {
		    max_dpix_dphi=dpix_dphi;
		}
		if(dddisp==0 && ddhdiv==0 && ddvdiv==0){
		    # movements due to "tangent mosaic" alone
		    if(dpix>max_dpix_tanm) max_dpix_tanm = dpix;
		    if(dphi>max_dphi_tanm) max_dphi_tanm = dphi;
		    if(dpix>max_pix_split){
			max_pix_split=dpix;
			pix_split_param="tanm";
		    }
		    if(dphi>max_phi_split){
			max_phi_split=dphi;
			phi_split_param="tanm";
		    }
		}
		if(ddtanm==0 && ddhdiv==0 && ddvdiv==0){
		    # movements due to dispersion alone
		    if(dpix>max_dpix_disp) max_dpix_disp = dpix;
		    if(dphi>max_dphi_disp) max_dphi_disp = dphi;
		    if(dpix>max_pix_split){
			max_pix_split=dpix;
			pix_split_param="disp";
		    }
		    if(dphi>max_phi_split){
			max_phi_split=dphi;
			phi_split_param="disp";
		    }
		}
		if(ddtanm==0 && dddisp==0 && ddvdiv==0){
		    # movements due to horiz divergecne alone
		    if(dpix>max_dpix_hdiv) max_dpix_hdiv = dpix;
		    if(dphi>max_dphi_hdiv) max_dphi_hdiv = dphi;
		    if(dpix>max_pix_split){
			max_pix_split=dpix;
			pix_split_param="hdiv";
		    }
		    if(dphi>max_phi_split){
			max_phi_split=dphi;
			phi_split_param="hdiv";
		    }
		}
		if(ddtanm==0 && dddisp==0 && ddhdiv==0){
		    # movements due to vertical divergecne alone
		    if(dpix>max_dpix_vdiv) max_dpix_vdiv = dpix;
		    if(dphi>max_dphi_vdiv) max_dphi_vdiv = dphi;
		    if(dpix>max_pix_split){
			max_pix_split=dpix;
			pix_split_param="vdiv";
		    }
		    if(dphi>max_phi_split){
			max_phi_split=dphi;
			phi_split_param="vdiv";
		    }
		}
	    }
	}
    }

    ++spot_corner;

}
}
}
}
#print "GOTHERE",h,k,l,1./d_star,max_dpix,max_dphi;
#print "GOTHERE",max_dpix_tanm,max_dpix_disp,max_dpix_hdiv,max_dpix_vdiv;
#print "GOTHERE",max_dphi_tanm,max_dphi_disp,max_dphi_hdiv,max_dphi_vdiv;
#print "GOTHERE",pix_split_param,phi_split_param;

dpix_kernel = speedup_pix;
dphi_kernel = speedup_phi*mosaic_spread/2+1e-30;

tan_steps_pix = int(max_dpix_tanm/dpix_kernel);
tan_steps_phi = int(max_dphi_tanm/dphi_kernel);
tan_steps = tan_steps_pix;
#if(tan_steps_phi > tan_steps) tan_steps = tan_steps_phi;
#phi_steps = tan_steps_phi;

phi_steps = int((max_dphi+mosaic_spread)/dphi_kernel);
#phi_steps = int((max_dpix_dphi*max_dphi)/dpix_kernel);
#print "GOTHERE", tan_steps_phi,tan_steps_pix;

disp_steps_pix = int(max_dpix_disp/dpix_kernel);
disp_steps_phi = int(max_dphi_disp/dphi_kernel);
disp_steps = disp_steps_pix;
if(disp_steps_phi > disp_steps) disp_steps = disp_steps_phi;

hdiv_steps_pix = int(max_dpix_hdiv/dpix_kernel);
hdiv_steps_phi = int(max_dphi_hdiv/dphi_kernel);
hdiv_steps = hdiv_steps_pix;
if(hdiv_steps_phi > hdiv_steps) hdiv_steps = hdiv_steps_phi;

vdiv_steps_pix = int(max_dpix_vdiv/dpix_kernel);
vdiv_steps_phi = int(max_dphi_vdiv/dphi_kernel);
vdiv_steps = vdiv_steps_pix;
if(vdiv_steps_phi > vdiv_steps) vdiv_steps = vdiv_steps_phi;



# there must be at least one sub-spot
if(phi_steps<1) phi_steps = 1
if(tan_steps<1) tan_steps = 1
if(disp_steps<1) disp_steps = 1
if(hdiv_steps<1) hdiv_steps = 1
if(vdiv_steps<1) vdiv_steps = 1
if(xtal_steps<1) xtal_steps = 1


# sanity check?
if(phi_steps>100) phi_steps = 100
if(tan_steps>100) tan_steps = 100
if(disp_steps>100) disp_steps = 100
if(hdiv_steps>100) hdiv_steps = 100
if(vdiv_steps>100) vdiv_steps = 100
if(xtal_steps>100) xtal_steps = 100

# for debugging
if(phimosaic_divisions) phi_steps =  phimosaic_divisions
if(tanmosaic_divisions) tan_steps =  tanmosaic_divisions
if(hdiv_divisions) hdiv_steps = hdiv_divisions
if(vdiv_divisions) vdiv_steps = vdiv_divisions
if(disp_divisions) disp_steps = disp_divisions

# pick "kernel"-sized chunks of each parameter
tanmosaic_step  = mosaic_spread/tan_steps;
phimosaic_step  = mosaic_spread/phi_steps;
dispersion_step = dispersion/disp_steps;
Hdiv_step       = Hdiv/hdiv_steps;
Vdiv_step       = Vdiv/vdiv_steps;
#xtal_step       = xtal_size/xtal_steps;

# inflate subrange size a bit to help blur them together?
sub_tanmosaic  = tanmosaic_step*1.0;
sub_phimosaic  = phimosaic_step*1.0;
sub_dispersion = dispersion_step*1.0;
sub_Hdiv       = Hdiv_step*1.0;
sub_Vdiv       = Vdiv_step*1.0;
#sub_xtal       = xtal_step*1.0;


#GOTHERE
#sub_phimosaic = 1e-6
#sub_tanmosaic = tanmosaic_step = 1e-6
#tanmosaic_step  = 0.6/RTD/tan_steps;
#phimosaic_step  = 0.2/RTD/phi_steps;


# iterate over range of mosaic spread and count up sub-rotations
mosaic_domains = 0;
mosweight_norm = 0;
for(phistep=1;phistep<=phi_steps;++phistep) {
for(tanstep=1;tanstep<=tan_steps;++tanstep) {
    phimosaic_dev  = phimosaic_step*(phistep-((phi_steps+1)/2))
    tanmosaic_dev  = tanmosaic_step*(tanstep-((tan_steps+1)/2))
    mosrad=sqrt(phimosaic_dev^2 + tanmosaic_dev^2)
    # constrain mosaic spread to a round spherical cap
    if(mosrad <= mosaic_spread/2){
	++mosaic_domains;
	# assume mosaic distribution obeys a Gaussian?
	mosweight_norm+=exp(-log(16)*(mosrad/(mosaic_spread+1e-30)*2)^2);
    }
}
}
if(mosweight_norm == 0) mosweight_norm = 1;


subspots = mosaic_domains*hdiv_steps*vdiv_steps*disp_steps

#print "GOTHERE", h,k,l,subspots,phi_steps,hdiv_steps,vdiv_steps,disp_steps,mosaic_domains, phi_range*RTD

# do this for speed?
if(subspots>max_subspots) {
#    print "GOTHERE",h,k,l,"is too big"
    continue;
}

# start generating sub-beams
for(dispstep=1;dispstep<=disp_steps;++dispstep) {
for(hstep=1;hstep<=hdiv_steps;++hstep) {
for(vstep=1;vstep<=vdiv_steps;++vstep) {
# and sub-crystals
for(phistep=1;phistep<=phi_steps;++phistep) {
for(tanstep=1;tanstep<=tan_steps;++tanstep) {

phimosaic_dev  = phimosaic_step*(phistep-((phi_steps+1)/2))
tanmosaic_dev  = tanmosaic_step*(tanstep-((tan_steps+1)/2))

# constrain mosaic spread to a round spherical cap
mosrad=sqrt(phimosaic_dev^2 + tanmosaic_dev^2)
if(mosrad > mosaic_spread/2) continue
mosweight = exp(-log(16)*(mosrad/(mosaic_spread+1e-30)*2)^2)/mosweight_norm;

# equal steps, but leave half a subdivision on either end
hdev       = Hdiv_step*(hstep-((hdiv_steps+1)/2))
vdev       = Vdiv_step*(vstep-((vdiv_steps+1)/2))
lambda_dev = dispersion_step*(dispstep-((disp_steps+1)/2))


# calculate the sub-spot position
XYphi(h,k,l,phimosaic_dev,tanmosaic_dev*mosaic_tanmult,lambda_dev,hdev,vdev);
if(Xdet=="nan") continue

# cannot hit the detector from behind!
if(cos(oblique_angle) <= 0) continue


if(det_size != "") {
    # don't bother if the ray will miss the detector completely
    if(Xdet < -10) continue
    if(Ydet < -10) continue
    if(Xdet*pixel > det_size+10) continue
    if(Ydet*pixel > det_size+10) continue
}

# compute reflection range (international tables C 2.2.7.3 p. 40)
sub_mosaic = sqrt(sub_tanmosaic^2 + sub_phimosaic^2)
# epsilon is supposed to be the radius of the relp, but our mosaicity is full width
epsilon = sin(twotheta)*( sub_phimosaic + sub_dispersion*sin(theta)/cos(theta))
phi_range = sqrt(((correlated_dispersion*d_star^2+zeta*sub_Hdiv)/Lorentz)^2 + sub_Vdiv^2) + epsilon/Lorentz


#phi_range = 1e-6
#print "GOTHERE", (max_dphi)*RTD, phi_range*RTD, sub_mosaic*RTD, mosaic_spread*RTD, mosaic_domains, 1./d_star;


# compute first-order approximation to spot shape



# angular spread of spot due to dispersion (in the radial direction)
del_theta_disp = sub_dispersion*sin(theta)/cos(theta)

# mosaic spread limits the radial spot size if dispersion is large?
#if(del_theta_disp > sub_tanmosaic) {
#    del_theta = sub_tanmosaic
#}else{
#    del_theta = del_theta_disp
#}
del_2theta = 2*del_theta_disp


# radial beam width at detector due to mosaic or dispersion (in mm)
mos_rwidth = airpath*sin(del_2theta)
# tangential beam width at detector due to mosaic only (in mm)
mos_twidth = airpath*(sub_tanmosaic)*2*sin(theta)



# broadening due to parallax at the sample
# beam width at detector due to illuminated path through crystal (in mm)
xtal_rwidth = beam_size*cos(twotheta)+xtal_size*sin(twotheta)
if(xtal_size<beam_size) {
    xtal_twidth = xtal_size
}else{
    xtal_twidth = beam_size
}

# total radial/tangential contribution to spot shape
Rsize = sqrt(mos_rwidth^2 + xtal_rwidth^2)/cos(oblique_angle)
Tsize = sqrt(mos_twidth^2 + xtal_twidth^2)


# spot size at detector due to divergence (in mm)
Hsize = airpath*sub_Hdiv/cos(Ydet_angle)
Vsize = airpath*sub_Vdiv/cos(Xdet_angle)


# calculate yield and effective phosphor depth
yield=1
impact(oblique_angle)
if(yield<=0) continue

# radial spot dimension will increase with parallax
apparent_Rsize = Rsize + det_thick_eff/1000*tan(oblique_angle)

# paralax will also shift the apparent radial spot position
paralax = 0.5*(apparent_Rsize - Rsize*cos(oblique_angle))

Xdet += paralax*radial_x;
Ydet += paralax*radial_y;

# correct for unwanted shift in beam center induced by spindle and detecctor tilts
#Xdet += dist/pixel*sin(detector_tilt)
#Ydet += dist/pixel*sin(spindle_swing+detector_twist)

# overall sub-spot shape will be the convolution of gaussians in these two directions
convolute_gauss(apparent_Rsize,Tsize,tangent_angle,Hsize,Vsize,0)
major_axis = wr
minor_axis = wt
final_angle = tilt


polar_correction();

#print "GOTHERE",phimosaic_dev,tanmosaic_dev,lambda_dev,hdev,vdev

#if(Xdet < 1500 || Ydet < 1500 || Xdet > 2000 || Ydet > 2000) continue

#if(! (h==3 && k==-2 && l==-2)) continue

printf "%4d %4d %4d   %11.4f %13.10f     %11.8g %3d    %11.4f %11.4f    %11.4f %11.4f %11.4f\n",\
          h,k,l, RTD*phi,RTD*phi_range, Lorentz*polar/yield/transmittance,subspots,  Xdet,Ydet,\
      RTD*final_angle,major_axis/pixel,minor_axis/pixel


}
}
}
}
}

}
}
}
}

}


################################################################################# 
#
#	model the interaction of an incident diffracted ray
# 	with a protective window and a detecting medium (phosphor) of finite depth
#
#	calculate the fraction of the maximum possible signal that is actually measured (yield)
#	as well as the "effective thickness" of the phosphor
#
############################################################################# 

function impact(beta) {

# effective attenuation depth is shallower at oblique incidence
mu_eff = muen_det*cos(beta)

# total energy deposited into the phosphor
yield = (1-safexp(-det_thick/mu_eff))*safexp(-window_thick/cos(beta)/mu_window)*safexp(-airpath/mu_air)

# however, some of the visible photons generated in the phosphor will get lost
# and this loss will increase with phosphor thickness
# especially if the mean-free-path of a visible photon in the phosphor (phosphor_mfp) is short
# but, this effect will also decrease with x-ray penetration depth into the phosphor
# the solution to this integral can be found in G. Lutz "Semiconductor Radiation Detectors" 1999
yield2 = (\
  mu_eff*phosphor_mfp/(mu_eff^2-phosphor_mfp^2)*(\
    (1+safexp(-det_thick/mu_eff))*coth(det_thick/phosphor_mfp)\
   -(1+safexp(-det_thick/mu_eff))/sinh(det_thick/phosphor_mfp)\
   -phosphor_mfp/mu_eff*(1-safexp(-det_thick/mu_eff))\
  )\
) * safexp(-window_thick/cos(beta)/mu_window)*safexp(-airpath/mu_air)

# use the more sophisticated model?
#yield = yield2

# the effective electro-optical (E-O) gain,
# provided we have a calibrated E-O gain at some wavelength
#EOgain = calibrated_EOgain*calbrated_EOgain_lambda/lambda*yield

# What depth would give us half the signal?

# if we use the not-so-fancy formula:
# 2*(1-safexp(-det_thick_eff/cos(beta)/muen_det)) = 1-safexp(-det_thick/cos(beta)/muen_det)
# det_thick_eff = -muen_det*cos(beta)*log((1+safexp(-det_thick*cos(beta)/muen_det))/2)

det_thick_eff = -mu_eff*log((1+safexp(-det_thick/mu_eff))/2)

}



################################################################################# 
#
#	misc useful maths...
#
############################################################################# 

function abs(x) {
  return sqrt(x*x)
}

function tan(angle) {
    return sin(angle)/cos(angle)
}

function safexp(x) { 
   if(x>700) return 1e300;
   if(x<-700) return 0;
   return exp(x);
}

function sinh(x) {
    return (safexp(x)-safexp(-x))/2;
}
function cosh(x) {
    return (safexp(x)+safexp(-x))/2;
}
function tanh(x) {
    return sinh(x)/cosh(-x);
}
function coth(x) {
    absx = abs(x)
    if(absx>350) return x/absx; 
    return 1./tanh(x);
}

# rotate a vector v_x,v_y,v_z about z then y then x axes
# creates new_x,new_y,new_z
function rotate(v_x,v_y,v_z,phix,phiy,phiz) {

new_x=v_x;new_y=v_y;new_z=v_z;


# rotate around z axis
rxx= cos(phiz); rxy=-sin(phiz); rxz= 0; 
ryx= sin(phiz); ryy= cos(phiz); ryz= 0;
rzx= 0;         rzy= 0;         rzz= 1;

rotated_x = new_x*rxx + new_y*rxy + new_z*rxz
rotated_y = new_x*ryx + new_y*ryy + new_z*ryz
rotated_z = new_x*rzx + new_y*rzy + new_z*rzz
new_x = rotated_x; new_y = rotated_y; new_z = rotated_z; 


# rotate around y axis
rxx= cos(phiy); rxy= 0;         rxz= sin(phiy); 
ryx= 0;         ryy= 1;         ryz= 0;
rzx=-sin(phiy); rzy= 0;         rzz= cos(phiy);

rotated_x = new_x*rxx + new_y*rxy + new_z*rxz
rotated_y = new_x*ryx + new_y*ryy + new_z*ryz
rotated_z = new_x*rzx + new_y*rzy + new_z*rzz
new_x = rotated_x; new_y = rotated_y; new_z = rotated_z; 


# rotate around x axis
rxx= 1;         rxy= 0;         rxz= 0; 
ryx= 0;         ryy= cos(phix); ryz=-sin(phix);
rzx= 0;         rzy= sin(phix); rzz= cos(phix);

rotated_x = new_x*rxx + new_y*rxy + new_z*rxz
rotated_y = new_x*ryx + new_y*ryy + new_z*ryz
rotated_z = new_x*rzx + new_y*rzy + new_z*rzz
new_x = rotated_x; new_y = rotated_y; new_z = rotated_z; 


}

################################################################################# 
#
#	convolution of two arbitrary Gaussians
#
#	returns: wr,wt,tilt
#
############################################################################# 
function convolute_gauss(wr1,wt1,tilt1,wr2,wt2,tilt2) {

  # convert input Gaussians into the general Gaussian equation:
  # exp(-(a^x^2 + b*x*y + c*y^2))
  # converting to reciprocal space, so we use 1/w

  a1 = (cos(  tilt1)*wr1)^2+(sin(  tilt1)*wt1)^2
  b1 = -sin(2*tilt1)*wr1^2  +sin(2*tilt1)*wt1^2
  c1 = (sin(  tilt1)*wr1)^2+(cos(  tilt1)*wt1)^2

  a2 = (cos(  tilt2)*wr2)^2+(sin(  tilt2)*wt2)^2
  b2 = -sin(2*tilt2)*wr2^2  +sin(2*tilt2)*wt2^2
  c2 = (sin(  tilt2)*wr2)^2+(cos(  tilt2)*wt2)^2

  # the convolution is a a product in reciprocal space, so we sum the parameters in the exponent
  a_conv = a1 + a2
  b_conv = b1 + b2
  c_conv = c1 + c2

  # decompose the general Gaussian equation parameters back into wr,wt and tilt
  #wr   = sqrt(0.5*(a_conv+c_conv - sqrt( a_conv^2 + b_conv^2 - 2*a_conv*c_conv + c_conv^2 )))
  #wt   = sqrt(0.5*(a_conv+c_conv + sqrt( a_conv^2 + b_conv^2 - 2*a_conv*c_conv + c_conv^2 )))
  #tilt = 0.5*atan2(-b_conv,a_conv-c_conv)
  # but we are going to check for singularities...
  dw=a_conv-c_conv
  if(dw==0 && b_conv == 0) {
   # difference in widths is zero, so it does not matter what the tilt is
   tilt = 0;
  } else {
   tilt = 0.5*atan2(-b_conv,dw);
  }
  square = a_conv^2 + b_conv^2 - 2*a_conv*c_conv + c_conv^2
  if(square < 0) {
    # how the fuck does this happen?
    print "GOTHERE1", wr1,wt1,tilt1,wr2,wt2,tilt2,"   ", square
  } else {
    square2 = 0.5*(a_conv+c_conv + sqrt( square ))
    if(square2 > 0) {
	wr = sqrt(square2)
    } else {
	# how the fick does this happen?
	print "GOTHERE2", wr1,wt1,tilt1,wr2,wt2,tilt2,"   ", square2
    }
    square3 = 0.5*(a_conv+c_conv - sqrt( square ))
    if(square3 > 0) {
	wt = sqrt(square3)
    } else {
	# how the fuck does this happen?
	print "GOTHERE3", wr1,wt1,tilt1,wr2,wt2,tilt2,"   ", square2
    }
  }

  if(tilt < -PI) tilt += 2*PI
  if(tilt >  PI) tilt -= 2*PI

#print "GOTHERE",wr1,wt1,tilt1*RTD,"   ",wr2,wt2,tilt2*RTD,"   ",wr,wt,tilt*RTD

}











function XYphi(h,k,l,phimos_dev,tanmos_dev,lambda_dev,hbeam_dev,vbeam_dev) {
############################################################################### 
#
#	calculate the intersection of the current relp
#	(defined by h,k,l and "sign")
#	with the Ewald sphere as it moves about the "phi" axis
#
#	modifies:
#	Xdet Ydet
#
#
#	variables used:
#	the A matrix (AX-CZ)
#	lambda_star0
#	spindle_swing
#
#
############################################################################### 

# default
if(sign+0==0) sign=1;

Xdet=Ydet=phi="nan"

# find the point in reciprocal space where this HKL resides at phi=0
x0=h*AX + k*BX + l*CX;
y0=h*AY + k*BY + l*CY;
z0=h*AZ + k*BZ + l*CZ;

# reciprocal-space magnitude of this HKL vector
d_star = sqrt(x0^2 + y0^2 + z0^2)


# apply the change of wavelength
lambda_star = lambda_star0 / (1 + lambda_dev)
lambda_star_sqr = lambda_star*lambda_star
lambda = 1/lambda_star

# apply the spindle swing
hbeam_dev = hbeam_dev + spindle_swing

# scattering length
#s = 2*d_star

# skip stuff outside the limiting sphere
if(d_star > 2*lambda_star) {
    # definitely outside Ewald sphere
    Xdet=Ydet="nan"

    # what could be done to bring this relp onto the Ewald sphere?
    required_lambda = 2/d_star;
    required_lambda_dev = required_lambda*lambda_star0 - 1;
    # there is no rotation that will bring this relp onto the Ewald sphere
    required_angle = "nan";

    return "outside";
}


# Bragg angle
sintheta = d_star/(2*lambda_star)
costheta = sqrt(1-sintheta^2)
theta    = atan2(sintheta,costheta)
twotheta = 2*theta



# apply the mosaic rotation along the phi rotation direction...

# rotate by phimos_dev around z axis
x1 = x0*cos(phimos_dev)-y0*sin(phimos_dev)
y1 = x0*sin(phimos_dev)+y0*cos(phimos_dev)


# apply the mosaic rotation normal to the phi rotation vector...

phi_radius0_sqr = (x0^2 + y0^2)
if(phi_radius0_sqr == 0) {
    # darn spherical cats...
    x0=1e-10
    phi_radius0_sqr = 1e-20
}
phi_radius0 = sqrt(phi_radius0_sqr)

# rotate the relp "toward" the z axis
# this move will never change tau
phi_radius = phi_radius0*cos(tanmos_dev) - z0*sin(tanmos_dev) 
phi_radius_sqr = phi_radius*phi_radius
x          = x1*phi_radius/phi_radius0
y          = y1*phi_radius/phi_radius0
z          = phi_radius0*sin(tanmos_dev) + z0*cos(tanmos_dev)


# forget the above ... use PURE rotations instead?

# asimuthal angle of the central relp (0 at z-pole)
alpha0 = atan2(phi_radius0,z0)
# "phi" rotation of central relp when "phi motor" = 0
tau0 = atan2(y0,x0)

# rotate a dummy relp (at the z-pole) by the prescribed mosaic rotations
rotate(0,0,d_star,phimos_dev,tanmos_dev,0);
x=rotated_x;
y=rotated_y;
z=rotated_z;

# now lay this down onto the location of the real relp
rotate(x,y,z,0,alpha0,0);
x=rotated_x;
y=rotated_y;
z=rotated_z;
rotate(x,y,z,0,0,tau0);
x=rotated_x;
y=rotated_y;
z=rotated_z;


phi_radius_sqr = (x^2 + y^2)
phi_radius = sqrt(phi_radius_sqr)


# find the values of phi that will put this HKL on the Ewald sphere...

# z is not going to change with phi


# calculate "phi" rotation of this HKL when "phi motor" = 0
tau = atan2(y,x)

# calculate the normalized (lambda=1) radius of the phi rotation circle for this hkl: xi (scribble)

# normalized (lambda=1) cyllindrical coordinates
xi_sqr = phi_radius_sqr/lambda_star_sqr
xi = sqrt(xi_sqr)
#xi = sqrt(x^2 + y^2)/lambda_star

# normalized coordinate of relp along rotation axis: zeta (squiggle)
zeta = z/lambda_star



# in case this relp never intersects the Ewald sphere ...
# what could be done to bring the phi circle of this relp to just touch the Ewald sphere?
# the declination of the relp from the z axis is alpha=atan2(xi,zeta)
# this vector must rotate away from the z axis until zeta=sin(twotheta)
# this vector must rotate away from the z axis until (90-theta)+hbeam_dev+alpha=90
# so, -required_angle = alpha-theta+hbeam_dev
if(xi!=0 && zeta!=0) {
    alpha=atan2(xi,zeta)
}else{
    alpha=0;
}
# required CHANGE in input mosaic or hdiv deviation to bring this hkl into the
# single point-of-contact condition
required_angle  = alpha-theta+hbeam_dev


# what wavelength will just touch the phi circle?
# 90-theta = 90-atan2(xi,zeta)
# theta = atan2(xi,zeta)
# lambda = 2*d*sin(atan2(xi,zeta)+hbeam_dev)
required_lambda = lambda
required_lambda_dev = "nan"
if(d_star != 0) {
    required_lambda = 2*sin(atan2(xi,zeta)+hbeam_dev)/d_star;
    required_lambda_dev = required_lambda*lambda_star0 - 1;
} else {
    required_lambda = lambda
    required_lambda_dev = "nan"
}


# where is the center of the Ewald sphere?  
# For normal beam geometry, it is at -lambda_star,0,0
# or xi,zeta,tau = 1,0,180

# what about not-quite-normal beam direction?
# Do "horiz" beam direction deviation (about "y") "first"

# quick rotate of (-1,0,0) around y axis
EC_x = -cos(hbeam_dev)
EC_y = 0
EC_z = sin(hbeam_dev)



# what is the radius of the (normalized) circle that is the intersection between the 
# (unit) Ewald sphere and the "z" level of this relp (zeta)?
# this will be independent of vdiv, since this is about the same axis as "phi"
# if hdev is positive, then EC_z is positive
# if hdev is large enough to make EC_z = zeta, then the center of the Ewald sphere
# is resting on the z-plane, and this radius is unity
r_Ez_sqr = 1 - (zeta - EC_z)^2



if(r_Ez_sqr <= 0) {
    # Ewald sphere never intersects this z-plane
    Xdet=Ydet="nan";
    return "outside";
}


# now, looking down the "z"-axis,
# the center of the unit Ewald sphere is at: -cos(hbeam_dev), 0
# where does the rotation circle of this hkl cross the Ewald sphere?
# call this intersection: xiP_x,xiP_y,zeta
# r_Ez^2 = (xiP_x+cos(hbeam_dev))^2 + xiP_y^2
# r_Ez^2 = xiP_x^2 +2*xiP_x*cos(hbeam_dev) +cos(hbeam_dev)^2 + xiP_y^2
# r_Ez^2 = xi^2 +2*xiP_x*cos(hbeam_dev) +cos(hbeam_dev)^2
# r_Ez^2 - xi^2 - cos(hbeam_dev)^2 = 2*xiP_x*cos(hbeam_dev)
# (r_Ez^2 - xi^2)/2/cos(hbeam_dev) - cos(hbeam_dev)/2 = xiP_x
xiP_x_novdev = (r_Ez_sqr - xi_sqr)/(2*cos(hbeam_dev)) - cos(hbeam_dev)/2


# solve for the two y coordinates (depends on "sign")
xiP_x_novdev_sqr = xiP_x_novdev*xiP_x_novdev
xiP_y_novdev_sqr = xi_sqr - xiP_x_novdev_sqr
if(xiP_y_novdev_sqr < 0) {
    # never intersects because xi is too small
    Xdet=Ydet="nan";

    # what "tanmosaic" rotation WOULD make these circles intersect?
    # this requires that the Ez circle rest on the "xi" circle, so: xiP_x=xi and xiP_y = 0
    # r_Ez^2 = (xiP_x+cos(hbeam_dev))^2 + xiP_y^2
    # r_Ez = xi+cos(hbeam_dev)
    # sqrt(1 - (zeta - sin(hbeam_dev))^2) = xi+cos(hbeam_dev)
    # hbeam_dev = Acos(r_Ez-xi)
#    cos_required_angle = sqrt(r_Ez_sqr)-xi
#    required_anglex = atan2(sqrt(1-cos_required_angle^2),cos_required_angle)
#print  "GOTHERE blind",xiP_y_novdev_sqr,cos_required_angle,required_anglex;

    # what energy would make these circles intersect?
    # r_Ez^2 = (xiP_x+cos(hbeam_dev))^2 + xiP_y^2
    # r_Ez = xi+cos(hbeam_dev)

    return "blind";
}
xiP_y_novdev     = sign*sqrt(xiP_y_novdev_sqr)



# now rotate the relp by vdev about the z axis

# rotate by vbeam_dev around z axis
xiP_x = xiP_x_novdev*cos(vbeam_dev)-xiP_y_novdev*sin(vbeam_dev)
xiP_y = xiP_x_novdev*sin(vbeam_dev)+xiP_y_novdev*cos(vbeam_dev)


# also rotate the Ewald sphere center
# quick result for rotation of hbeam_dev about "y" then vbeam_dev about "z" axis
EC_x = -cos(hbeam_dev)*cos(vbeam_dev)
EC_y = -cos(hbeam_dev)*sin(vbeam_dev)
EC_z = sin(hbeam_dev)



# calculate absolute rotation of this intersection point of this HKL with the Ewald sphere
intercept_rot = atan2(xiP_y,xiP_x)

# actual "phi" rotation that will bring x,y onto xp,yp will depend on
# the phi=0 setting (tau) and the "vertical" beam direction (rotation about the spindle)
phi = (intercept_rot - tau)
if(phi < -PI) phi += 2*PI
if(phi >  PI) phi -= 2*PI


# what is the unit vector representing the diffracted beam?


# the diffracted ray simply connects the center of the Ewald sphere to the relp
diffray_x = xiP_x - EC_x
diffray_y = xiP_y - EC_y
diffray_z = zeta - EC_z





# Lorentz correction of normal beam rotation camera geometry
#Lorentz_sqr0 = sin(twotheta)^2 - (zeta)^2
#Lorentz0 = 0
#if(Lorentz_sqr0>0) Lorentz0 = sqrt(Lorentz_sqr0)

# the Lorentz correction is the component of the velocity of the spot (in lambda=1 units)
# along the diffracted beam path
# the Lorentz factor L is the inverse of the Lorentz correction

# what is the vector representing the velocity of the relp?
# there is an identical triangle to xp,yp,phi_radius that is rotated 90 degrees about the relp
# here the phi_radius becomes tangent to the circle of rotation and equal to the velocity
velocity_x = xiP_y
velocity_y = -xiP_x
velocity_z = 0

#speed = sqrt(velocity_x^2 + velocity_y^2 + velocity_z^2)
speed = xi

# the Lorentz correction is the dot product of the unit diffracted ray and this velocity
Lorentz = diffray_x*velocity_x + diffray_y*velocity_y + diffray_z*velocity_z
if(Lorentz==0){
    # wow, what are the chances of that!
    #return "grazing"
}

# absolute value
Lorentz_sqr = Lorentz*Lorentz
Lorentz = sqrt(Lorentz_sqr)


# how fast will a small phi move change the spot position?
sqr = speed^2 - Lorentz_sqr
if(sqr>0){
    tan_speed = sqrt(sqr);
}else{
    tan_speed = 0;
}


# evaluate the absorber envelope
min_exitpath = 1e99;
min_entrypath = 1e99;
mu_exitpath = 1e99;
mu_entrypath = 1e99;
for(i=1;i<=absorbers;++i) {
    # rotate this absorber facet to the phi value where this spot diffracts
    rotate(absorber_x[i],absorber_y[i],absorber_z[i],0,0,phi);
    # find the component of the absorber vector in the direction of the diffracted ray
    proj=new_x*diffray_x + new_y*diffray_y + new_z*diffray_z;

    # negative projection means the facet is more than 90 degrees away
    if(proj>0) {

        # ratio of the pathlength through material to the slab depth is the same as the
        # ratio of the slab depth to the projection of the depth vector onto the beam path
        pathlength = (absorber_x[i]*absorber_x[i]+absorber_y[i]*absorber_y[i]+absorber_z[i]*absorber_z[i])/proj;
    }else{
	pathlength = 1e99;
    }
    
    # the shortest path is the path out of the convex solid
    if(pathlength<min_exitpath) {
	min_exitpath = pathlength;
	mu_exitpath = absorber_mu[i];
    }

    # find the component of the absorber vector in the direction of the incident ray
    proj=new_x*EC_x + new_y*EC_y + new_z*EC_z;

    # negative projection means the facet is more than 90 degrees away
    if(proj>0) {

        # ratio of the pathlength through material to the slab depth is the same as the
        # ratio of the slab depth to the projection of the depth vector onto the beam path
        pathlength = (absorber_x[i]*absorber_x[i]+absorber_y[i]*absorber_y[i]+absorber_z[i]*absorber_z[i])/proj;
    }else{
	pathlength = 1e99;
    }
    
    # the shortest path is the path out of the convex solid
    if(pathlength<min_entrypath) {
	min_entrypath = pathlength;
	mu_entrypath = absorber_mu[i];
    }
}
if(min_entrypath == 1e99) min_entrypath=0;
if(min_exitpath == 1e99) min_exitpath=0;
exponent = min_entrypath/mu_entrypath+min_exitpath/mu_exitpath;
if(exponent<100){
    transmittance=exp(-exponent);
}else{
    # path is effectively infinite?
    transmittance=1e-99;
}



# where will this spot hit the detector?


# transform the diffracted ray into the detector coordinate system
xd = detector_X_x*diffray_x + detector_X_y*diffray_y + detector_X_z*diffray_z
yd = detector_Y_x*diffray_x + detector_Y_y*diffray_y + detector_Y_z*diffray_z
zd = detector_Z_x*diffray_x + detector_Z_y*diffray_y + detector_Z_z*diffray_z

# where does the central direct-beam (1,0,0 rotated by spindle_swing) hit?
db_x = cos(spindle_swing)
db_y = 0
db_z = -sin(spindle_swing)

xd0 = detector_X_x*db_x + detector_X_y*db_y + detector_X_z*db_z
yd0 = detector_Y_x*db_x + detector_Y_y*db_y + detector_Y_z*db_z
zd0 = detector_Z_x*db_x + detector_Z_y*db_y + detector_Z_z*db_z

# shortest distance from crystal to detector plane
#xtf=dist*zd0

# convert to pixel coordinates
Xdet = (xd/zd*dist -xd0/zd0*dist +xbeam )/pixel + 1
Ydet = (yd/zd*dist -yd0/zd0*dist +ybeam )/pixel

dXdet_dphi = tan_speed*dist/pixel

# make a memorable unit vector in the radial direciton
radial_x = xd
radial_y = yd


# how far will this beam travel through the air?
airpath = dist/zd

# fractional radial distance from direct beam on detector surface
rd = sqrt(xd*xd + yd*yd)

# obliquity: what is the angle of incidence of this ray with the detector surface?
oblique_angle = atan2(rd,zd)
# component in the "x" direction
Xdet_angle    = atan2(xd,zd)
# component in the "y" direction
Ydet_angle    = atan2(yd,zd)

# what is the angle with the detector X-axis?
tangent_angle = atan2(xd,yd)

return 1;
}




function polar_correction(gravity_x,gravity_y,gravity_z) {
#################################################################
#
#	compute the polarization correction
#
# polarization correction requires: psi
# this is the angle between the "vertical" axis (normal to E-vector and beam) 
# and the projection of the vector normal to the reflecting plane onto the plane normal to the beam

# the normal to the reflecting plane should always be contained within the plane containing the incident and
# reflected beams, so the projection of the vector normal to the reflecting plane onto the plane normal to
# the beam should overlap the projection of the diffracted ray onto the plane normal to the incident beam

# assume that the polarization is always horizontal, and it is the spindle that is messed up?

# unit vector in the direction of the incident beam
incident_x = -EC_x
incident_y = -EC_y
incident_z = -EC_z

# arbitrary unit vector
if(gravity_x == "") {
    gravity_x = 0
    gravity_y = 1
    gravity_z = 0
}

# the E-vector must always be normal to the direction of propagation
Evector_x = incident_y*gravity_z - incident_z*gravity_y
Evector_y = incident_z*gravity_x - incident_x*gravity_z
Evector_z = incident_x*gravity_y - incident_y*gravity_x

# "vertical" axis that is normal to the incident beam and the E-vector
Bvector_x = Evector_y*incident_z - Evector_z*incident_y
Bvector_y = Evector_z*incident_x - Evector_x*incident_z
Bvector_z = Evector_x*incident_y - Evector_y*incident_x

# get components of diffracted ray projected onto the E-B plane
diffray_Evector = diffray_x*Evector_x + diffray_y*Evector_y + diffray_z*Evector_z
diffray_Bvector = diffray_x*Bvector_x + diffray_y*Bvector_y + diffray_z*Bvector_z

# compute the polarization angle
psi = -atan2(diffray_Bvector,diffray_Evector)
twopsi = 2*psi

# correction for polarized incident beam
polar = 2/(1 + cos(twotheta)^2 - polar_fac*cos(2*psi)*sin(twotheta)^2)

return polar;
}
