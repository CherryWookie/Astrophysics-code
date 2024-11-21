# Hello world this is my code

import numpy as np
import pylab as pl
import astropy
import galpy
import astropy.units as u
from galpy.orbit import Orbit
import galpy.util.coords as coords
from galpy.potential import MWPotential2014

N = 1000    # Number of objects

# Cartesian reference system (kpc):
X = 8
Y = 0
Z = 0.025

# DISK galaxies:
# Mass Distribution: p(R,Z) = p(Rsolar,0)e^(Rsolar/L) * e^(-R/L -(Z + Zsolar)/H)
# Ligh profile: I(R) = I_0 * e^-(R/h_R)


# Physical Paramters
R0 = 8*u.kpc #kpc
z0 = 0.025*u.kpc #kpc
V0 = 220*u.km/u.s   #km/s
u_sun = -11.1 #km/s
v_sun = 7.25*u.km / u.s  #km/s
w_sun = 12.24 #km/s


np.random.seed(0)

# for local files, just use this function #
R = np.genfromtxt('Gaia_100pcWD.csv', delimiter=',', names=True, dtype=None, encoding="ASCII")
print(R.dtype)

# masking the galactic components by using the lable defined by Torres et al.
mask_thin = np.where(R['Ilabel_0thin1thick2halo']==0)
mask_thick = np.where(R['Ilabel_0thin1thick2halo']==1)
mask_halo = np.where(R['Ilabel_0thin1thick2halo']==2)



star_thin = R[mask_thin]
star_thick = R[mask_thick]
star_halo = R[mask_halo]


ts = np.linspace(0,0.125,10000)*u.Gyr

o_thin = Orbit([star_thin['RA_deg'][0],star_thin['DEC_deg'][0],1/star_thin['parallax_mas'][0],
                star_thin['pmramasyr'][0],star_thin['pmdec_masyr'][0], 0.],
                 radec=True, ro=R0, zo=z0, vo=V0, solarmotion=v_sun)
o_thin.integrate(ts,MWPotential2014)

o_thick = Orbit([star_thick['RA_deg'][10],star_thick['DEC_deg'][10],1/star_thick['parallax_mas'][10],
                star_thick['pmramasyr'][10],star_thick['pmdec_masyr'][10], 0.],
                 radec=True, ro=R0, zo=z0, vo=V0, solarmotion=v_sun)
o_thick.integrate(ts,MWPotential2014)

o_halo = Orbit([star_halo['RA_deg'][-1],star_halo['DEC_deg'][-1],1/star_halo['parallax_mas'][-1],
                star_halo['pmramasyr'][-1],star_halo['pmdec_masyr'][-1], 0.],
                 radec=True, ro=R0, zo=z0, vo=V0, solarmotion=v_sun)
o_halo.integrate(ts,MWPotential2014)

### integrate all ###
a_thin = Orbit([star_thin['RA_deg'][0:5000],star_thin['DEC_deg'][0:5000],1/star_thin['parallax_mas'][0:5000],
                star_thin['pmramasyr'][0:5000],star_thin['pmdec_masyr'][0:5000], np.zeros(len(star_thin[0:5000]))],
                 radec=True, ro=R0, zo=z0, vo=V0, solarmotion=v_sun)
a_thin.integrate(ts,MWPotential2014)

a_thick = Orbit([star_thick['RA_deg'],star_thick['DEC_deg'],1/star_thick['parallax_mas'],
                star_thick['pmramasyr'],star_thick['pmdec_masyr'], np.zeros(len(star_thick))],
                 radec=True, ro=R0, zo=z0, vo=V0, solarmotion=v_sun)
a_thick.integrate(ts,MWPotential2014)

a_halo = Orbit([star_halo['RA_deg'],star_halo['DEC_deg'],1/star_halo['parallax_mas'],
                star_halo['pmramasyr'],star_halo['pmdec_masyr'], np.zeros(len(star_halo))],
                 radec=True, ro=R0, zo=z0, vo=V0, solarmotion=v_sun)
a_halo.integrate(ts,MWPotential2014)


# printing eccentricity and angular momentum (z component)
print('e = ', round(o_thin.e(),3), 'Lz = ', o_thin.Lz(quantity=True))
print('e = ', round(o_thick.e(),3), 'Lz = ', o_thick.Lz(quantity=True))
print('e = ', round(o_halo.e(),3), 'Lz = ', o_halo.Lz(quantity=True))


fig = pl.figure(figsize=(7,7))
fig.subplots_adjust(wspace=0.3)
ax1 = fig.add_subplot(111)

ax1.plot(star_halo['V_kms'],np.sqrt(star_halo['U_kms']**2 + star_halo['W_kms']), marker='.', ms=3, ls='none', label='halo')
ax1.plot(star_thick['V_kms'],np.sqrt(star_thick['U_kms']**2 + star_thick['W_kms']), marker='.', ms=3, ls='none', label='thick disk')
ax1.plot(star_thin['V_kms'],np.sqrt(star_thin['U_kms']**2 + star_thin['W_kms']), marker='.', ms=3, ls='none', label='thin disk')

ax1.set_xlabel(r'$V_{\rm LSR}$ [km/s]')
ax1.set_ylabel(r'$\sqrt{U^2_{\rm LSR} + W^2_{\rm LSR}}$ [km/s]')

pl.legend(frameon=False)


pl.show()