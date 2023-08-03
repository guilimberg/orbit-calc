import agama
import math
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
from pygaia.astrometry.coordinates import CoordinateTransformation, Transformations
from pygaia.astrometry.vectorastrometry import cartesian_to_spherical, astrometry_to_phase_space,phase_space_to_astrometry

import time

start = time.time()
end = time.time()
#print(end - start)

#Set the units and the potential
agama.setUnits(mass=1,length=1,velocity=1)
potential = agama.Potential("McMillan17.ini")

def ICRS_to_GAL(phi, theta, muphistar, mutheta):
	"""
	phi       - The longitude-like angle of the position of the source (radians)
	theta     - The latitude-like angle of the position of the source (radians)
	muphistar - Value of the proper motion in the longitude-like angle, multiplied by cos(latitude).
	mutheta   - Value of the proper motion in the latitude-like angle.
	"""
	# transformation ICRS to GAL
	ctICRS2GAL =CoordinateTransformation(Transformations.ICRS2GAL)
    	#
	ell,         bee = ctICRS2GAL.transform_sky_coordinates(phi,theta)
	muellstar, mubee = ctICRS2GAL.transform_proper_motions (phi,theta, muphistar,mutheta)

	return ell,bee,muellstar,mubee

def astrometric_2_cartesian(DM, ell_radian, bee_radian, HRV, muellstar, mubee):
	parallax_mas = np.power(10., 2.-DM/5.)
	xHelio_pc, yHelio_pc, zHelio_pc, vxHelio, vyHelio, vzHelio = astrometry_to_phase_space(ell_radian, bee_radian, parallax_mas, muellstar, mubee, HRV)
	xHelio_kpc = xHelio_pc/1000.
	yHelio_kpc = yHelio_pc/1000.
	zHelio_kpc = zHelio_pc/1000.
	return xHelio_kpc, yHelio_kpc, zHelio_kpc, vxHelio, vyHelio, vzHelio

def cartesian_2_cylindrical(x,y,z,vx,vy,vz):
	R =np.sqrt(x**2 + y**2)
	r =np.sqrt(x**2 + y**2 + z**2)

	phi=np.arctan2(y,x)
	vR     = (x*vx + y*vy)/R
	vTHETA = (R*vR - z*vz)/r
	vPHI   = (x*vy - y*vx)/R

	return R, phi, vR, vTHETA, vPHI

#calculate orbital parameters
def calculate_orbital_parameters(xGC_norm, yGC_norm, zGC_norm, vxGC_norm, vyGC_norm, vzGC_norm, inttime):
    	# galpy convention
	x_agama  =   xGC_norm
	y_agama  =   yGC_norm
	z_agama  =   zGC_norm
	vx_agama =   vxGC_norm
	vy_agama =   vyGC_norm
	vz_agama =   vzGC_norm
	
	R_agama   =  np.sqrt(x_agama**2 + y_agama**2)
	vR_agama  =  ( x_agama*vx_agama + y_agama*vy_agama )/R_agama
	vT_agama  = -( x_agama*vy_agama - y_agama*vx_agama )/R_agama
	phi_agama =  np.arctan2(y_agama,x_agama)


	inttime=20./0.98
	numsteps=5000
	times = np.linspace(0, inttime, numsteps)
	times_c, c_orb_car = agama.orbit(ic=[x_agama,y_agama,z_agama,vx_agama,vy_agama,vz_agama], potential=potential, time=inttime, trajsize=numsteps)

	x = c_orb_car[:,0]
	y = c_orb_car[:,1]
	z = c_orb_car[:,2]
	vx= c_orb_car[:,3]
	vy= c_orb_car[:,4]
	vz= c_orb_car[:,5]
	#Lx= np.array(c_orb_car[:,1])*np.array(c_orb_car[:,5]) - np.array(c_orb_car[:,2])*np.array(c_orb_car[:,4])
	#Ly= np.array(c_orb_car[:,2])*np.array(c_orb_car[:,3]) - np.array(c_orb_car[:,0])*np.array(c_orb_car[:,5])
	#Lz= np.array(c_orb_car[:,0])*np.array(c_orb_car[:,4]) - np.array(c_orb_car[:,1])*np.array(c_orb_car[:,3])	
	R = np.sqrt(x**2 + y**2)
	r = np.sqrt(x**2 + y**2 + z**2)

	rmin = np.min(r)
	rmax = np.max(r)
	Rmax = np.max(R)
	zmax = np.max(np.fabs(z))
	ecc  = (rmax-rmin)/(rmax+rmin)

	#Lx   = Lx[0] 
	#Ly   = Ly[0]
	#Lz   = Lz[0]

	#vx   = vx[0] 
	#vy   = vy[0]
	#vz   = vz[0]

	#xgal = x[0] 
	#ygal = y[0]
	#zgal = z[0]

	return rmin, rmax, Rmax, zmax, ecc

"""
#position and velocity of the Sun ## BOVY 2020 ## based on Gaia EDR3
_xGC_sun_kpc = -8.224
_yGC_sun_kpc = 0.
_zGC_sun_kpc = 0.
_vxGC_sun    = 11.10
_vyGC_sun    =  7.20 + 243. # Vsun + Vlsr
_vzGC_sun    =  7.25
"""

# position and velocity of the Sun
# this parameters were chosen to match Anna's email
_xGC_sun_kpc = -8.2 # McMillan2017 + Bland-HawthornGehrard2016 
_yGC_sun_kpc =  0.
_zGC_sun_kpc =  0. # ????
_vxGC_sun    = 11.10  # Sun's peculiar motion from Schronrich2011
_vyGC_sun    = 12.24 + 232.8 # Vsun + Vlsr (lsr from mcmillan)
_vzGC_sun    =  7.25


#Open the data file
data = []
with open("your_data_here",'r') as file:
	data = file.readlines()
#Obtain the RA, DEC, distance, radial velocity, pmRA, and pmDEC from the file
name = []
ra = []
dec = []
dist = []
rv = []
pmra = []
pmdec = []
e_pmra = []
e_pmdec = []
e_dist = []
e_HRV = [] 
counter = -1
for line in data:
	counter = counter + 1
	if counter != 0:
		split_line = line.split(',')
		name.append(split_line[0])				#1
		ra.append(split_line[1])				#2
		dec.append(split_line[2])				#3
		dist.append(float(split_line[3]))		#4
		rv.append(float(split_line[4]))			#5
		pmra.append(float(split_line[5]))		#6
		pmdec.append(float(split_line[6]))		#7
		e_pmra.append(float(split_line[7]))		#8
		e_pmdec.append(float(split_line[8]))	#9
		e_dist.append(float(split_line[9]))		#10
		e_HRV.append(float(split_line[10]))		#11

#Put the coordinates into a bin to get them easier in degrees
#coordinates = SkyCoord(ra_dec,unit=(u.hourangle,u.deg))
print("Os dados já foram carregados! Que alegria! :D")
for i in range (len(ra)):
	ra[i] = float(ra[i])

for i in range (len(dec)):
	dec[i] = float(dec[i])

#Put the stellar position and velocities into an object
icrs_coord = coord.ICRS(ra=ra*u.deg, 
			dec=dec*u.deg,
		        distance=dist*u.kpc,
			pm_ra_cosdec=pmra*u.mas/u.yr,
			pm_dec=pmdec*u.mas/u.yr,
			radial_velocity=rv*u.km/u.s)

#Convert to galactic coordinates
galac_coord = icrs_coord.transform_to(coord.Galactocentric)

#Make an array with the cartesian positions and velocities of all the stars
galac_points = np.zeros(shape=(len(galac_coord.x),6))
for i in np.arange(0,len(galac_coord.x)):
	galac_points[i][0] = galac_coord.x.value[i]
	galac_points[i][1] = galac_coord.y.value[i]
	galac_points[i][2] = galac_coord.z.value[i]
	galac_points[i][3] = galac_coord.v_x.value[i]
	galac_points[i][4] = galac_coord.v_y.value[i]
	galac_points[i][5] = galac_coord.v_z.value[i]

#Calculate the potential and kinetic energy
ke = 0.5*np.sum(galac_points[:,3:]**2, axis=1)
pe = potential.potential(galac_points[:,0:3])
galac_energy = ke + pe

#Find the action using the potential
act_finder = agama.ActionFinder(potential)
galac_action = act_finder(galac_points)

j_r = []
j_z = []
j_phi = []
for action in galac_action:
	j_r.append(action[0])
	j_z.append(action[1])
	j_phi.append(action[2])

x = []
y = []
v_x = []
v_y = []
v_z = []
for star in galac_points:
	x.append(star[0])
	y.append(star[1])
	v_x.append(star[3])
	v_y.append(star[4])
	v_z.append(star[5])

j_x = []
j_y = []
v_r = []
v_phi = []
for i in np.arange(0,len(v_x)):
	R = np.sqrt(x[i]**2+y[i]**2)
	j_x.append(j_r[i]*np.cos(j_phi[i]))
	j_y.append(j_r[i]*np.sin(j_phi[i]))
	v_r.append((x[i]*v_x[i]+y[i]*v_y[i])/R)
	v_phi.append((x[i]*v_y[i]-y[i]*v_x[i])/R)

#Convert the degrees to radians for the position of the stars and gather data
#into convenient lists
RA_rad = []
DEC_rad = []
HRV_kms = []
pmRA = []
pmDEC = []
parallax_mas = []
count = -1
for star in icrs_coord:
	count += 1
	RA_rad.append(star.ra.deg*np.pi/180.0)
	DEC_rad.append(star.dec.deg*np.pi/180.0)
	HRV_kms.append(star.radial_velocity.value)
	pmRA.append(star.pm_ra_cosdec.value)
	pmDEC.append(star.pm_dec.value)
	parallax_mas.append(1.0/dist[count])

#Find the error for the proper motion based on Gaia DR2 error (Took G < 15 off the website)
#Find the error for the distance based on Gaia DR2 error (Took off the website)
"""
e_pmRAstar_masyr = np.ones(len(name))
#e_dist = np.ones(len(name))
for i in np.arange(0,len(name)):
	G = Gmag[i]
	if (G < 15.0):
		e_pmRAstar_masyr[i] = 0.06
		#e_dist[i] = (0.04)**2/parallax_mas[i]
	elif (15.0 < G <=17.0):
		e_pmRAstar_masyr[i] = 0.06 + (G-15.)*(0.2-0.06)/2.
		#e_dist[i] = (0.1)**2/parallax_mas[i]
	elif (17. < G):
		e_pmRAstar_masyr[i] = 0.20 + (G-17.)*(1.2-0.2)/3.
		#e_dist[i] = (0.7)**2/parallax_mas[i]
"""
e_pmRAstar_masyr = e_pmra
e_pmDECstar_masyr = e_pmdec
###############################################################################
#Calculate the actions and velocities and their errors using monte carlo methods

#Set the length and of the simulation
num_stars = len(name)
_N_MonteCarlo = 100#1000
np.random.seed(666)

#Arrays to put the results of the montecarlo simulation into
JrJzJphiE_Nstar_Nmontecarlo = np.ones((4, num_stars, _N_MonteCarlo))
vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo = np.ones((20, num_stars, _N_MonteCarlo))

Jr_median = []
Jr_std = []
Jr_16 = []
Jr_84 = []

Jz_median = []
Jz_16 = []
Jz_84 = []

Jphi_median = []
Jphi_16 = []
Jphi_84 = []

Jx_mean = []
Jx_std = []

Jy_mean = []
Jy_std = []

vr_median = []
vr_std = []
vr_16 = []
vr_84 = []

vx_mean = []
vx_std = []

vy_mean = []
vy_std = []

vz_median = []
vz_16 = []
vz_84 = []

vphi_median = []
vphi_16 = []
vphi_84 = []

energy_median = []
energy_16 = []
energy_84 = []

phi_mean = []
phi_std = []

rperi_median = []
rperi_16 = []
rperi_84 = []

rapo_median = []
rapo_16 = []
rapo_84 = []

Rmax_median = []
Rmax_16 = []
Rmax_84 = []

zmax_median = []
zmax_16 = []
zmax_84 = []

ecc_median = []
ecc_mean = []
ecc_16 = []
ecc_84 = []

Lx_median = []
Lx_16 = []
Lx_84 = []

Ly_median = []
Ly_16 = []
Ly_84 = []

Lz_median = []
Lz_16 = []
Lz_84 = []

Ltot_median = []
Ltot_16 = []
Ltot_84 = []

Lperp_median = []
Lperp_16 = []
Lperp_84 = []

incl_median = []
incl_16 = []
incl_84 = []

xgal_median = []
xgal_16 = []
xgal_84 = []
ygal_median = []
ygal_16 = []
ygal_84 = []
zgal_median = []
zgal_16 = []
zgal_84 = []

vx_median = []
vx_16 = []
vx_84 = []
vy_median = []
vy_16 = []
vy_84 = []


for i in np.arange(0,num_stars):
	end = time.time()
	print("Estrelinha:", i, "| Tempo:", round((end - start)/3600,2), "hora(s)" )
	#print(i)
	for j in np.arange(0,_N_MonteCarlo):
		RA_rad_ij = RA_rad[i] #no error
		DEC_rad_ij = DEC_rad[i] #no error
		HRV_kms_ij = HRV_kms[i] + e_HRV[i] * np.random.normal() #put in random error
		parallax_mas_i = parallax_mas[i]
		mean_pmRA_pmDEC = [pmRA[i],pmDEC[i]]
		#Get the error for RA and DEC
		e_pmRA = e_pmRAstar_masyr[i]
		e_pmDEC = e_pmDECstar_masyr[i]
		#Calculate the sigma from RA and DEC
		Sigma_pmRA = e_pmRA**2
		Sigma_pmRA_pmDEC = e_pmRA*e_pmDEC #Upper limit
		Sigma_pmDEC = e_pmDEC**2
		Sigma = [[Sigma_pmRA, Sigma_pmRA_pmDEC],
			[Sigma_pmRA_pmDEC, Sigma_pmDEC]]
		#Create an array of the RA and DEC
		pmRAstar_masyr_ij_array, pmDECstar_masyr_ij_array = np.random.multivariate_normal(mean=mean_pmRA_pmDEC, cov=Sigma, size=1).T
		pmRAstar_masyr_ij = pmRAstar_masyr_ij_array[0]
		pmDECstar_masyr_ij = pmDECstar_masyr_ij_array[0]
		#Create a distance with error
		dist_kpc_ij = dist[i] + np.random.normal()*e_dist[i]
		while dist_kpc_ij < 0.0:
			dist_kpc_ij = dist[i] + np.random.normal()*e_dist[i]
		parallax_mas_ij = 1.0/dist_kpc_ij
		#Calculate a magnitude from the calculated parallax
		DM_mag_ij = 5.*(np.log10(1./parallax_mas_ij) + 2.0)
		#Put the star into galactic coordinates
		l_rad_ij, b_rad_ij, pmlstar_masyr_ij, pmbstar_masyr_ij = ICRS_to_GAL(RA_rad_ij, DEC_rad_ij, pmRAstar_masyr_ij, pmDECstar_masyr_ij)
		#Convert to Cartesian
		xH, yH, zH, vxH, vyH, vzH = astrometric_2_cartesian(DM_mag_ij, l_rad_ij, b_rad_ij, HRV_kms_ij, pmlstar_masyr_ij, pmbstar_masyr_ij)
		#Normalize the position
		xGC_norm = (xH + _xGC_sun_kpc)
		yGC_norm = (yH + _yGC_sun_kpc)
		zGC_norm = (zH + _zGC_sun_kpc)
 		#normalize  velocity
		vxGC_norm = (vxH + _vxGC_sun)
		vyGC_norm = (vyH + _vyGC_sun)
		vzGC_norm = (vzH + _vzGC_sun)
 		#calculate angular momentum vector components
		LxGC_norm = yGC_norm*vzGC_norm - zGC_norm*vyGC_norm
		LyGC_norm = zGC_norm*vxGC_norm - xGC_norm*vzGC_norm
		LzGC_norm = xGC_norm*vyGC_norm - yGC_norm*vxGC_norm
		#calculate total angular momentum, orthogonal component of angular momentum, and orbital inclination
		Ltot = (LxGC_norm**2 + LyGC_norm**2 + LzGC_norm**2)**0.5
		Lperp = (LxGC_norm**2 + LyGC_norm**2)**0.5
		incl = (np.arccos(LzGC_norm/Ltot)) * (180/(np.pi))
		#cylindrical coordinate
		R, phi, vR, vTHETA, vPHI = cartesian_2_cylindrical(xGC_norm, yGC_norm, zGC_norm, vxGC_norm, vyGC_norm, vzGC_norm)
		vPERP = np.sqrt(vR**2 + vzGC_norm**2)
		#Calculate the energy
		energy = potential.potential([xGC_norm,yGC_norm,zGC_norm]) + 0.5*(vxGC_norm**2+vyGC_norm**2+vzGC_norm**2)
		#Calculate the actions
		xv6d = np.array( [xGC_norm, yGC_norm, zGC_norm, vxGC_norm, vyGC_norm, vzGC_norm]).T
		Jr, Jz, Jphi = act_finder(xv6d)

		rmin, rmax, Rmax, zmax, ecc = -1., -1., -1., -1., -1.
		if (np.isnan(Jr)):
			Jr = -1.
		if (np.isnan(Jz)):
			Jz = -1.
		if (np.isnan(energy)):
			print(energy, parallax_mas_ij)
			sys.exit(1)
		#I'm pretty sure I do not need this right now, but maybe later
		#Get the orbital parameters
		if (energy<0.):
			rmin,rmax,Rmax,zmax,ecc = calculate_orbital_parameters(xGC_norm,yGC_norm,zGC_norm,vxGC_norm,vyGC_norm,vzGC_norm,10)
			#print(i, j)
		#else:
			#print("here")
		#Put the actions and the energy into an array
		JrJzJphiE_Nstar_Nmontecarlo[0,i,j] = Jr
		JrJzJphiE_Nstar_Nmontecarlo[1,i,j] = Jz
		JrJzJphiE_Nstar_Nmontecarlo[2,i,j] = Jphi
		JrJzJphiE_Nstar_Nmontecarlo[3,i,j] = energy
		#Put the velocities into an array
		#I put phi into VPERP because I didn't want to change all the VPERP
		vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[0,i,j] = vR
		vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[1,i,j] = vzGC_norm
		vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[2,i,j] = vPHI
		vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[3,i,j] = phi
		#Put the distances and eccentricities into an array
		vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[4,i,j] = rmin
		vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[5,i,j] = rmax
		vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[6,i,j] = Rmax
		vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[7,i,j] = zmax
		vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[8,i,j] = ecc
		vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[9,i,j] = LxGC_norm
		vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[10,i,j] = LyGC_norm
		vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[11,i,j]= LzGC_norm
		vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[12,i,j]= Ltot
		vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[13,i,j]= Lperp
		vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[14,i,j]= incl
		vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[15,i,j]= xGC_norm
		vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[16,i,j]= yGC_norm
		vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[17,i,j]= zGC_norm
		vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[18,i,j]= vxGC_norm
		vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[19,i,j]= vyGC_norm

	#Calculate the median and std of the actions
	Jr_median.append(np.median(JrJzJphiE_Nstar_Nmontecarlo[0,i]))
	Jr_std.append(np.std(JrJzJphiE_Nstar_Nmontecarlo[0,i]))
	Jr_16.append(np.percentile(JrJzJphiE_Nstar_Nmontecarlo[0,i], 16))
	Jr_84.append(np.percentile(JrJzJphiE_Nstar_Nmontecarlo[0,i], 84))
	Jz_median.append(np.median(JrJzJphiE_Nstar_Nmontecarlo[1,i]))
	Jz_16.append(np.percentile(JrJzJphiE_Nstar_Nmontecarlo[1,i], 16))
	Jz_84.append(np.percentile(JrJzJphiE_Nstar_Nmontecarlo[1,i], 84))
	Jphi_median.append(np.median(JrJzJphiE_Nstar_Nmontecarlo[2,i]))
	Jphi_16.append(np.percentile(JrJzJphiE_Nstar_Nmontecarlo[2,i], 16))
	Jphi_84.append(np.percentile(JrJzJphiE_Nstar_Nmontecarlo[2,i], 84))
        #Calculate the median and std of the energy
	energy_median.append(np.median(JrJzJphiE_Nstar_Nmontecarlo[3,i]))
	energy_16.append(np.percentile(JrJzJphiE_Nstar_Nmontecarlo[3,i], 16))
	energy_84.append(np.percentile(JrJzJphiE_Nstar_Nmontecarlo[3,i], 84))
	#Calculate the median and std of the velocities
	vr_median.append(np.median(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[0,i]))
	vr_std.append(np.std(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[0,i]))
	vr_16.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[0,i], 16))
	vr_84.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[0,i], 84))
	vz_median.append(np.median(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[1,i]))
	vz_16.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[1,i], 16))
	vz_84.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[1,i], 84))
	vphi_median.append(np.median(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[2,i]))
	vphi_16.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[2,i], 16))
	vphi_84.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[2,i], 84))
	rperi_median.append(np.median(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[4,i]))
	rperi_16.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[4,i], 16))
	rperi_84.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[4,i], 84))
	rapo_median.append(np.median(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[5,i]))
	rapo_16.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[5,i], 16))
	rapo_84.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[5,i], 84))
	Rmax_median.append(np.median(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[6,i]))
	Rmax_16.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[6,i], 16))
	Rmax_84.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[6,i], 84))
	zmax_median.append(np.median(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[7,i]))
	zmax_16.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[7,i], 16))
	zmax_84.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[7,i], 84))
	ecc_median.append(np.median(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[8,i]))
	ecc_16.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[8,i], 16))
	ecc_84.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[8,i], 84))
	ecc_mean.append(np.mean(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[8,i]))

	Lx_median.append(np.median(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[9,i]))
	Lx_16.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[9,i], 16))
	Lx_84.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[9,i], 84))

	Ly_median.append(np.median(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[10,i]))
	Ly_16.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[10,i], 16))
	Ly_84.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[10,i], 84))

	Lz_median.append(np.median(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[11,i]))
	Lz_16.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[11,i], 16))
	Lz_84.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[11,i], 84))
	
	Ltot_median.append(np.median(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[12,i]))
	Ltot_16.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[12,i], 16))
	Ltot_84.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[12,i], 84))
	
	Lperp_median.append(np.median(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[13,i]))
	Lperp_16.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[13,i], 16))
	Lperp_84.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[13,i], 84))

	incl_median.append(np.median(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[14,i]))
	incl_16.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[14,i], 16))
	incl_84.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[14,i], 84))

	xgal_median.append(np.median(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[15,i]))
	xgal_16.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[15,i], 16))
	xgal_84.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[15,i], 84))

	ygal_median.append(np.median(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[16,i]))
	ygal_16.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[16,i], 16))
	ygal_84.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[16,i], 84))

	zgal_median.append(np.median(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[17,i]))
	zgal_16.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[17,i], 16))
	zgal_84.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[17,i], 84))

	vx_median.append(np.median(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[18,i]))
	vx_16.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[18,i], 16))
	vx_84.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[18,i], 84))

	vy_median.append(np.median(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[19,i]))
	vy_16.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[19,i], 16))
	vy_84.append(np.percentile(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[19,i], 84))

	#Calculate the median and std of phi (Used to convert to Cartesian Coordinates)
	phi_mean.append(np.mean(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[3,i]))
	phi_std.append(np.std(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[3,i]))
	"""
	#Calculate the median and std of Cartesian actions and velocities
	Jx_mean.append(np.mean(Jr_median[i]*np.cos(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[3,i])))
	Jx_std.append(np.sqrt(np.cos(phi_mean[i])**2*Jr_std[i]*2+Jr_median[i]**2*np.sin(phi_mean[i])**2*phi_std[i]**2))
	Jy_mean.append(np.mean(Jr_median[i]*np.sin(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[3,i])))
	Jy_std.append(np.sqrt(np.sin(phi_mean[i])**2*Jr_std[i]**2+Jr_median[i]**2*np.cos(phi_mean[i])**2*phi_std[i]**2))


	vx_mean.append(np.mean(vr_median[i]*np.cos(vr_vz_vphi_vPERP_rPERI_rAPO_zMAX_ecc__Nstar__NMonteCarlo[3,i])))
	vx_std.append(np.sqrt(np.cos(phi_mean[i])**2*vr_std[i]*2+vr_median[i]**2*np.sin(phi_mean[i])**2*phi_std[i]**2))

	vy_mean.append(vr_median[i]*np.sin(-vphi_median[i]))
	vy_std.append(np.sqrt(np.sin(phi_mean[i])**2*vr_std[i]**2+vr_median[i]**2*np.cos(phi_mean[i])**2*phi_std[i]**2))
	"""

with open('output_orb_param.csv','w') as file:
	file.write('Name,v_r16,v_r50,v_r84,v_phi16,v_phi50,v_phi84,v_z16,v_z50,v_z84,j_r16,j_r50,j_r84,j_phi16,j_phi50,j_phi84,j_z16,j_z50,j_z84,energy16,energy50,energy84,rperi16,rperi50,rperi84,rapo16,rapo50,rapo84,Rmax16,Rmax50,Rmax84,zmax16,zmax50,zmax84,ecc16,ecc50,ecc84,ecc_mean,Lx16,Lx50,Lx84,Ly16,Ly50,Ly84,Lz16,Lz50,Lz84,Ltot16,Ltot50,Ltot84,Lperp16,Lperp50,Lperp84,incl16,incl50,incl84,xgal16,xgal50,xgal84,ygal16,ygal50,ygal84,zgal16,zgal50,zgal84,vx16,vx50,vx84,vy16,vy50,vy84,vz16,vz50,vz84\n')
	for i in np.arange(0,len(name)):
	
		#Write the name of the object
		file.write(name[i]+',')
		
		#Write the velocities in CYLINDRICAL fram
		file.write(str(vr_16[i])+',')
		file.write(str(vr_median[i])+',')		
		file.write(str(vr_84[i])+',')

		file.write(str(vphi_16[i])+',')
		file.write(str(vphi_median[i])+',')
		file.write(str(vphi_84[i])+',')

		file.write(str(vz_16[i])+',')
		file.write(str(vz_median[i])+',')
		file.write(str(vz_84[i])+',')
		
		#Write the ACTION-ANGLE variables
		file.write(str(Jr_16[i])+',')
		file.write(str(Jr_median[i])+',')
		file.write(str(Jr_84[i])+',')

		file.write(str(Jphi_16[i])+',')
		file.write(str(Jphi_median[i])+',')
		file.write(str(Jphi_84[i])+',')

		file.write(str(Jz_16[i])+',')
		file.write(str(Jz_median[i])+',')
		file.write(str(Jz_84[i])+',')

		#Write the ORBITAL ENERGY
		file.write(str(energy_16[i])+',')
		file.write(str(energy_median[i])+',')
		file.write(str(energy_84[i])+',')

		#Write the Apocentric/Pericentric Distance, z max and the eccentricity
		file.write(str(rperi_16[i])+',')
		file.write(str(rperi_median[i])+',')
		file.write(str(rperi_84[i])+',')

		file.write(str(rapo_16[i])+',')
		file.write(str(rapo_median[i])+',')
		file.write(str(rapo_84[i])+',')

		file.write(str(Rmax_16[i])+',')
		file.write(str(Rmax_median[i])+',')
		file.write(str(Rmax_84[i])+',')

		file.write(str(zmax_16[i])+',')
		file.write(str(zmax_median[i])+',')
		file.write(str(zmax_84[i])+',')

		file.write(str(ecc_16[i])+',')
		file.write(str(ecc_median[i])+',')
		file.write(str(ecc_84[i])+',')
		file.write(str(ecc_mean[i])+',')

		#Write the cartesian components of ANGULAR MOMENTUM and orbital INCLINATION
		file.write(str(Lx_16[i])+',')
		file.write(str(Lx_median[i])+',')
		file.write(str(Lx_84[i])+',')

		file.write(str(Ly_16[i])+',')
		file.write(str(Ly_median[i])+',')
		file.write(str(Ly_84[i])+',')

		file.write(str(Lz_16[i])+',')
		file.write(str(Lz_median[i])+',')
		file.write(str(Lz_84[i])+',')
		
		file.write(str(Ltot_16[i])+',')
		file.write(str(Ltot_median[i])+',')
		file.write(str(Ltot_84[i])+',')

		file.write(str(Lperp_16[i])+',')
		file.write(str(Lperp_median[i])+',')
		file.write(str(Lperp_84[i])+',')

		file.write(str(incl_16[i])+',')
		file.write(str(incl_median[i])+',')
		file.write(str(incl_84[i])+',')

		#Write the cartesian POSITIONS
		file.write(str(xgal_16[i])+',')
		file.write(str(xgal_median[i])+',')
		file.write(str(xgal_84[i])+',')

		file.write(str(ygal_16[i])+',')
		file.write(str(ygal_median[i])+',')
		file.write(str(ygal_84[i])+',')

		file.write(str(zgal_16[i])+',')
		file.write(str(zgal_median[i])+',')
		file.write(str(zgal_84[i])+',')

		#Write the cartesian VELOCITIES
		file.write(str(vx_16[i])+',')
		file.write(str(vx_median[i])+',')
		file.write(str(vx_84[i])+',')

		file.write(str(vy_16[i])+',')
		file.write(str(vy_median[i])+',')
		file.write(str(vy_84[i])+',')

		file.write(str(vz_16[i])+',')
		file.write(str(vz_median[i])+',')
		file.write(str(vz_84[i])+'\n')


