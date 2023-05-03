import numpy as np
import matplotlib.pyplot as plt

mixfrac0 = np.load('slices2d/00000_rhoMixFrac.npy',allow_pickle=True).item()
rhoh0 = np.load('slices2d/00000_rhoh.npy',allow_pickle=True).item()
density0 = np.load('slices2d/00000_density.npy',allow_pickle=True).item()

h_min, h_max = np.amin(rhoh0['data']/density0['data']),np.amax(rhoh0['data']/density0['data'])
mf_min, mf_max = np.amin(mixfrac0['data']/density0['data']),np.amax(mixfrac0['data']/density0['data'])

nt = 100

mixfrac = np.load(f'slices2d/{nt:05d}_rhoMixFrac.npy',allow_pickle=True).item()
rhoh = np.load(f'slices2d/{nt:05d}_rhoh.npy',allow_pickle=True).item()
density = np.load(f'slices2d/{nt:05d}_density.npy',allow_pickle=True).item()

mid = int(mixfrac['data'].shape[1]/2)

mf = (mixfrac['data'][:,mid]/density['data'][:,mid] - mf_min) / (mf_max-mf_min)
h =  (rhoh['data'][:,mid]/density['data'][:,mid] - h_min) / (h_max-h_min)

plt.plot(mixfrac['y'],mf,label='Z')
plt.plot(rhoh['y'],h,label='H',ls='--')
plt.legend()
plt.show()
