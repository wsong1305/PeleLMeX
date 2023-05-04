import numpy as np
import matplotlib.pyplot as plt

mixfrac0_0 = np.load('slices2d/00000_rhoMixFrac0.npy',allow_pickle=True).item()
rhoh_0 = np.load('slices2d/00000_rhoh.npy',allow_pickle=True).item()
density_0 = np.load('slices2d/00000_density.npy',allow_pickle=True).item()

h_min, h_max = np.amin(rhoh_0['data']/density_0['data']),np.amax(rhoh_0['data']/density_0['data'])
mf0_min, mf0_max = np.amin(mixfrac0_0['data']/density_0['data']),np.amax(mixfrac0_0['data']/density_0['data'])

nt = 80

mixfrac0 = np.load(f'slices2d/{nt:05d}_rhoMixFrac0.npy',allow_pickle=True).item()
mixfrac1 = np.load(f'slices2d/{nt:05d}_rhoMixFrac1.npy',allow_pickle=True).item()
rhoh = np.load(f'slices2d/{nt:05d}_rhoh.npy',allow_pickle=True).item()
density = np.load(f'slices2d/{nt:05d}_density.npy',allow_pickle=True).item()

mid = int(mixfrac0['data'].shape[1]/2)

mf0 = (mixfrac0['data'][:,mid]/density['data'][:,mid] - mf0_min) / (mf0_max-mf0_min)
mf1 = (mixfrac1['data'][:,mid]/density['data'][:,mid] - mf0_min) / (mf0_max-mf0_min)
h =  (rhoh['data'][:,mid]/density['data'][:,mid] - h_min) / (h_max-h_min)

plt.plot(mixfrac0['y'],mf0,label='Z_0')
plt.plot(mixfrac1['y'],mf1,label='Z_1')
plt.plot(rhoh['y'],h,label='H',ls='--')
plt.legend()
plt.show()
