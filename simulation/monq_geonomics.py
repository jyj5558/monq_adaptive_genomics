#ml anaconda
#ml use.own
#ml conda-env/geonomics-py3.11.9
#ml proj
#ml geos
#ml gdal
#unset LD_LIBRARY_PATH
#unset PYTHONPATH
#conda activate geonomics

#python

import geonomics as gnx
import matplotlib.pyplot as plt

# set Matplotlib's default plotting style and plot size
plt.rcParams["figure.figsize"] = (9,4)

# make our params file 
gnx.make_parameters_file('monq_params.py',
                         layers=[{'type':'file', 'change':True}]*1,
                         species=[{'movement':True, 'movement_surface':True,
                                   'genomes':True, 'n_traits':1}],
                         data=True, stats=True)

#  then use it to make a model
mod = gnx.make_model('./monq_gnx_params.py', verbose=True)

# save model object
import pickle

with open("monq_gnx_mod.pkl", "wb") as m:
  pickle.dump(mod, m)
  
# take a look at the resulting object
mod

# interactively plot the resulting object
mod.plot(lyr=0, spp=0)

# run the burn-in, then plot phenotypes
mod.walk(T=10000, mode='burn')
mod.plot(lyr=0, spp=0)

fig = plt.figure(figsize=(18,36))
ax1 = fig.add_subplot(111)
mod.plot_phenotype(spp=0, trt=0)
plt.show() 

# run the first 5k, then plot phenotypes
mod.walk(5000, 'main')
mod.plot(lyr=0, spp=0)

fig = plt.figure(figsize=(18,36))
ax1 = fig.add_subplot(111)
mod.plot_phenotype(spp=0, trt=0)
plt.show() 