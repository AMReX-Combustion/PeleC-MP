#!/usr/bin/env python
#
# Generate a table of the 2D velocity fluctuations for the homogeneous
# isotropic turbulence case for specific values of u' and L11
# NC - Grid size in wavespace
# NX - Grid size in physical space
#
# Order of operations:
#   1. velocity fluctuations generated on a NC^2 grid in wavenumber space
#   2. Coefficients associated to wavenumbers that cannot be represented on
#      the desired grid are set to 0 (sharp wavenumber cutoff)
#   3. inverse Fourier transform of the velocity fluctuations (NC^2 grid)
#   4. velocity fluctuations resampled on the desired grid (NX^2)
#

# ========================================================================
#
# Imports
#
# ========================================================================
import argparse
import sys
import time
from datetime import timedelta
import numpy as np
import scipy.interpolate as spi
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt


# ========================================================================
#
# Parse arguments
#
# ========================================================================
parser = argparse.ArgumentParser(
    description="Generate the velocity fluctuations for the HIT IC"
)
parser.add_argument(
    "-k0", help="Wave number containing highest energy", type=float, default=-1.
) 
parser.add_argument("-N", help="Resolution", type=int, default=16)
parser.add_argument(
    "-s", "--seed", help="Random number generator seed", type=int, default=42
)
parser.add_argument(
    "-up", help="Turbulent fluctations", type=float, default=1.
)
parser.add_argument(
    "-l11", help="Velocity autocorrelation lenght-scale", type=float, default=0.044
)
parser.add_argument(
    "-L", help="Physical domain length", type=float, default=1.
)
args = parser.parse_args()

# ===============================================================================
#
# Some defaults variables
#
# ===============================================================================
plt.rc("text", usetex=True)
plt.rc("font", family="serif", serif="Times")
cmap_med = [
    "#F15A60",
    "#7AC36A",
    "#5A9BD4",
    "#FAA75B",
    "#9E67AB",
    "#CE7058",
    "#D77FB4",
    "#737373",
]
cmap = [
    "#EE2E2F",
    "#008C48",
    "#185AA9",
    "#F47D23",
    "#662C91",
    "#A21D21",
    "#B43894",
    "#010202",
]
dashseq = [
    (None, None),
    [10, 5],
    [10, 4, 3, 4],
    [3, 3],
    [10, 4, 3, 4, 3, 4],
    [3, 3],
    [3, 3],
]
markertype = ["s", "d", "o", "p", "h"]

# ========================================================================
#
# Function definitions
#
# ========================================================================
def div0(a, b):
    """ Ignore division by 0, just replace it by 0,

    From: http://stackoverflow.com/questions/26248654/numpy-return-0-with-divide-by-zero
    e.g. div0( [-1, 0, 1], 0 ) -> [0, 0, 0]
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        c = np.true_divide(a, b)
        c[~np.isfinite(c)] = 0  # -inf inf NaN
    return c


# ========================================================================
def abs2(x):
    """This is equivalent to np.abs(x)**2 or x*np.conj(x)

    To make it faster, add this right before the function definition
    import numba
    @numba.vectorize([numba.float64(numba.complex128),numba.float32(numba.complex64)])
    """
    return x.real ** 2 + x.imag ** 2


# ========================================================================
#
# Main
#
# ========================================================================

# Timer
start = time.time()

# ========================================================================
# 1. velocity fluctuations generated on a NC^2 grid in wavenumber space

# Dimension of the large cube
N = 512
#N = 256
halfN = int(round(0.5 * N))
xs = 0
xe = args.L
dx = args.L / N
if (args.k0 < 0.):
    args.k0 = 8./3.*np.sqrt(2./np.pi)/args.l11

# Only work if N and args.N are even
if not ((args.N % 2 == 0) and N % 2 == 0):
    print("N or args.N is not even. Exiting")
    sys.exit(1)

amp = 2.*np.pi/args.L
# Get cell centered values and meshed grid
x = np.linspace(xs, xe, N + 1)
xc = (x[1:] + x[:-1]) / 2.  # get cell center coordinates
X, Y = np.meshgrid(xc, xc, indexing="ij")
# Get the wave numbers and associated quantities
k = np.concatenate((np.arange(halfN), np.arange(-halfN, 0, 1)), axis=0)
khalf = np.arange(halfN + 1)
k1, k2 = np.meshgrid(khalf, k, indexing="ij")
kmag = np.sqrt((amp*k1) ** 2 + (amp*k2) ** 2)

# Energy spectrum
Ek = (
    32./3.*np.sqrt(2./np.pi)*
    (kmag**4)/(args.k0**5)*
    np.exp(-2.*(kmag**2)/(args.k0**2))
)

# Draw random numbers
np.random.seed(args.seed)
phi = np.random.uniform(0, 2 * np.pi, np.shape(kmag))
eps_new = 2.*np.pi/args.L
eps = args.k0/10000.

# uft = amp*np.sqrt(div0(Ek*(amp*k2)**2, np.pi * (kmag ** 3)))*np.exp(1j * phi)
# vft = -div0((amp*k1)**2*uf, amp*k2)

ps = np.exp(1j * phi)
# Just to get the correct shape
uf = np.zeros_like(ps)
vf = np.zeros_like(ps)
# Initialize velocities in spectral space
for i in range(0,halfN):
    for j in range(0,N):
        kx = amp*k1[i,j]
        ky = amp*k2[i,j]
        if (abs(ky) > eps):
            uf[i,j] = amp*np.sqrt(Ek[i,j]*ky**2/(np.pi*kmag[i,j]**3))*ps[i,j]
            if (abs(kx) > eps):
                vf[i,j] = -kx**2*uf[i,j]/ky
        elif (abs(kx) > eps):
            vf[i,j] = amp*np.sqrt(Ek[i,j]*kx**2/(np.pi*kmag[i,j]**3))*ps[i,j]

# Impose the 3D spherical symmetry (to ensure we have a real signal)
# equiv: uf[-l,-m,0] = np.conj(uf[ l, m,0]) for l=0..N/2 and m=0..N/2
uf[0, N:halfN:-1] = np.conj(uf[0, 1:halfN])
vf[0, N:halfN:-1] = np.conj(vf[0, 1:halfN])

# Normalize. Because we are generating the data in wavenumber space,
# we have to multiply by N**3 because in the definition of the numpy
# ifftn there is a 1/N**n.
uf = uf * N ** 2
vf = vf * N ** 2

# # Quick check on energy content (make sure you add both the current
# # contribution and the one we are neglecting because we are assuming
# # real input data)
# print('Energy = int E(k) dk = 0.5 * int (uf**2 + vf**2 wf**2) dk1 dk2 dk3 = {0:.10f} ~= 3/2'.format(
#     (np.sum(abs2(uf          ) + abs2(vf          ) + abs2(wf          )) +
#      np.sum(abs2(uf[:,:,1:-1]) + abs2(vf[:,:,1:-1]) + abs2(wf[:,:,1:-1])))
#     * 0.5 / N**6))

# ========================================================================
# 2. Coefficients associated to wavenumbers that cannot be represented
# on the desired grid are set to 0 (sharp wavenumber cutoff)
kmagc = 0.5 * args.N
uf[kmag > kmagc] = 0.0
vf[kmag > kmagc] = 0.0

# ========================================================================
# 3. inverse Fourier transform of the velocity fluctuations (NC^2 grid)
u = np.fft.irfftn(uf, s=(N, N))
v = np.fft.irfftn(vf, s=(N, N))

# Another energy content check
print(
    "Energy = 1/V * int E(x,y) dV = 0.5/V * int (u**2 + v**2) dx dy = {0:.10f} ~= 1".format(
        np.sum(u ** 2 + v ** 2) * 0.5 * (dx / args.L) ** 2/args.up
    )
)

# ========================================================================
# 4. velocity fluctuations re-sampled on the desired grid (N^3)
xr = np.linspace(xs, xe, args.N + 1)
xrc = (xr[1:] + xr[:-1]) / 2.
Xr, Yr = np.meshgrid(xrc, xrc, indexing="ij")

Xr = Xr.reshape(-1, order="F")
Yr = Yr.reshape(-1, order="F")

ur = spi.interpn((xc, xc), u, (Xr, Yr), method="linear")*args.up
vr = spi.interpn((xc, xc), v, (Xr, Yr), method="linear")*args.up

# ========================================================================
# Save the data in Fortran ordering
fname = "hit_ic_{0:d}_{1:d}.dat".format(int(args.k0), args.N)
data = np.vstack((Xr, Yr, ur, vr)).T
np.savetxt(fname, data, fmt="%.18e", delimiter=",", header="x, y, u, v")

# output timer
end = time.time() - start
print("Elapsed time " + str(timedelta(seconds=end)) + " (or {0:f} seconds)".format(end))
