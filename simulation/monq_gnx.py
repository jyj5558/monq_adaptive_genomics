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
import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend
import matplotlib.pyplot as plt
import io
import sys

# make our params file 
#gnx.make_parameters_file('monq_params.py',
#                         layers=[{'type':'file', 'change':True}]*1,
#                         species=[{'movement':True, 'movement_surface':True,
#                                   'genomes':True, 'n_traits':1}],
#                         data=True, stats=True)

# then use it to make a model
mod = gnx.make_model('./monq_gnx_params.py', verbose=True)
#mod = gnx.make_model('./monq_test_params.py', verbose=True)

# take a look at the resulting object
mod

# interactively plot the resulting object
mod.plot(lyr=0, spp=0)
fig = plt.gcf()
fig.savefig("monq_mod_start.png", dpi=300, bbox_inches='tight')
plt.close(fig)

# run the burn-in, then plot phenotypes
mod.walk(T=10000, mode='burn')

species = list(mod.comm.values())[0]
nonneut_array = species.gen_arch.nonneut_loci 

# Unpacking the four values
nonneut_1, nonneut_2, nonneut_3, nonneut_4 = nonneut_array[:4]

print(nonneut_1)
print(nonneut_2)
print(nonneut_3)
print(nonneut_4)

random_individs = mod.get_random_individs(n=10000, spp=0)

# set Matplotlib's default plotting style and plot size
plt.rcParams["figure.figsize"] = (9,4)

mod.plot(lyr=0, spp=0, individs=random_individs, color='none', edge_color='black', size=0.1, alpha=0.5)
fig = plt.gcf()
fig.savefig("monq_mod_burnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_phenotype(individs=random_individs, edge_color='face', size=0.1, alpha=1)
fig = plt.gcf()
fig.savefig("monq_phenotype_burnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

#density raster is not being drawn
#mod.plot_density(spp=0, individs=random_individs, size=0.1)
#fig = plt.gcf()
#fig.savefig("monq_density_burnin.png", dpi=300, bbox_inches='tight')
#plt.close(fig)

mod.plot_genetic_PCA(individs=random_individs, edge_color='face', size=0.1, alpha=1)
fig = plt.gcf()
fig.savefig("monq_pca_burnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_allele_frequencies(spp=0, color='red')
fig = plt.gcf()
fig.savefig("monq_af_burnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_hist_fitness()
fig = plt.gcf()
fig.savefig("monq_hist_ft_burnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_movement_surface(style='vect')
fig = plt.gcf()
fig.savefig("monq_movement_sf_burnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_dispersal_surface(style='vect')
fig = plt.gcf()
fig.savefig("monq_dispersal_sf_burnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_dispersal(n_individs=10000, include_start_points=False, size=0.1)
fig = plt.gcf()
fig.savefig("monq_dispersal_burnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_movement(n_individs=10000, include_start_points=False, size=0.1)
fig = plt.gcf()
fig.savefig("monq_movement_burnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_demographic_pyramid()
fig = plt.gcf()
fig.savefig("monq_pyramid_burnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_pop_growth(expected=True, actual=True, expected_color='red', actual_color='blue')
fig = plt.gcf()
fig.savefig("monq_growth_burnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_genotype(locus=nonneut_1, individs=random_individs, edge_color='face', size=0.1, alpha=1) #non-neutral locus from non-neutral data csv file
fig = plt.gcf()
fig.savefig("monq_genotype" + str(nonneut_1) + "_burnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_genotype(locus=nonneut_2, individs=random_individs, edge_color='face', size=0.1, alpha=1) #non-neutral locus from non-neutral data csv file
fig = plt.gcf()
fig.savefig("monq_genotype" + str(nonneut_2) + "_burnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_genotype(locus=nonneut_3, individs=random_individs, edge_color='face', size=0.1, alpha=1) #non-neutral locus from non-neutral data csv file
fig = plt.gcf()
fig.savefig("monq_genotype" + str(nonneut_3) + "_burnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_genotype(locus=nonneut_4, individs=random_individs, edge_color='face', size=0.1, alpha=1) #non-neutral locus from non-neutral data csv file
fig = plt.gcf()
fig.savefig("monq_genotype" + str(nonneut_4) + "_burnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

plt.rcParams["figure.figsize"] = (6,6) 
mod.plot_fitness(individs=random_individs, edge_color='face', size=0.5, alpha=0.5, fit_cmap='seismic')
fig = plt.gcf()
plt.subplots_adjust(right=0.8) 
fig.savefig("monq_fitness_burnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

# run the preliminary
#mod.walk(100, 'main')

# figure out non-neutral burn-in
#import numpy as np
#count = 0
#while count <= 150:
#  mod.walk(1, 'main')
#  arr = mod.get_fitness()
#  new_arr = np.insert(arr, 0, count)
#  if count == 0:
#    np.savetxt("fitness_trace.csv", new_arr[np.newaxis, :], delimiter=",", fmt="%f")
#  else:
#    with open("fitness_trace.csv", "a") as f:
#      np.savetxt(f, new_arr[np.newaxis, :], delimiter=",", fmt="%f")
#  count += 1

#back to linux terminal  
#awk '{file=sprintf("fitness_%d.csv", int((NR-1)/10)); print > file}' fitness_trace.csv
#for j in {0..10}; do awk -F',' '{sum=0; count=0; for (i=2; i<=NF; i++) {if ($i != "") {sum+=$i; count++}} print sum/count}' fitness_$j.csv >> fit_avg.csv; done

#in R
#library(tseries)
#file_prefix <- "fit_avg_"
#file_extension <- ".csv"
#for (i in 0:10) {
#    # Construct file name
#    file_name <- paste0(file_prefix, i, file_extension)
    
#    # Read the single-column time series
#    ts_data <- scan(file_name, quiet = TRUE)  # scan() reads a single-column file directly
    
#    # Run ADF test
#    adf_result <- adf.test(ts_data, k = 10)
    
#    # Print results
#    cat("\nADF Test Result for", file_name, ":\n")
#    print(adf_result)
#}

# Non-neutral runs
count = 0
while count < 60:
  mod.walk(10, 'main')
  count += 10

random_individs = mod.get_random_individs(n=10000, spp=0)

# set Matplotlib's default plotting style and plot size
plt.rcParams["figure.figsize"] = (9,4)

mod.plot(lyr=0, spp=0, individs=random_individs, color='none', edge_color='black', size=0.1, alpha=0.5)
fig = plt.gcf()
fig.savefig("monq_mod_nntrburnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_phenotype(individs=random_individs, edge_color='face', size=0.1, alpha=1)
fig = plt.gcf()
fig.savefig("monq_phenotype_nntrburnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

#density raster is not being drawn
#mod.plot_density(spp=0, individs=random_individs, size=0.1)
#fig = plt.gcf()
#fig.savefig("monq_density_nntrburnin.png", dpi=300, bbox_inches='tight')
#plt.close(fig)

mod.plot_genetic_PCA(individs=random_individs, edge_color='face', size=0.1, alpha=1)
fig = plt.gcf()
fig.savefig("monq_pca_nntrburnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_allele_frequencies(spp=0, color='red')
fig = plt.gcf()
fig.savefig("monq_af_nntrburnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_hist_fitness()
fig = plt.gcf()
fig.savefig("monq_hist_ft_nntrburnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_movement_surface(style='vect')
fig = plt.gcf()
fig.savefig("monq_movement_sf_nntrburnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_dispersal_surface(style='vect')
fig = plt.gcf()
fig.savefig("monq_dispersal_sf_nntrburnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_dispersal(n_individs=10000, include_start_points=False, size=0.1)
fig = plt.gcf()
fig.savefig("monq_dispersal_nntrburnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_movement(n_individs=10000, include_start_points=False, size=0.1)
fig = plt.gcf()
fig.savefig("monq_movement_nntrburnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_demographic_pyramid()
fig = plt.gcf()
fig.savefig("monq_pyramid_nntrburnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_pop_growth(expected=True, actual=True, expected_color='red', actual_color='blue')
fig = plt.gcf()
fig.savefig("monq_growth_nntrburnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_genotype(locus=nonneut_1, individs=random_individs, edge_color='face', size=0.1, alpha=1) #non-neutral locus from non-neutral data csv file
fig = plt.gcf()
fig.savefig("monq_genotype" + str(nonneut_1) + "_nntrburnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_genotype(locus=nonneut_2, individs=random_individs, edge_color='face', size=0.1, alpha=1) #non-neutral locus from non-neutral data csv file
fig = plt.gcf()
fig.savefig("monq_genotype" + str(nonneut_2) + "_nntrburnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_genotype(locus=nonneut_3, individs=random_individs, edge_color='face', size=0.1, alpha=1) #non-neutral locus from non-neutral data csv file
fig = plt.gcf()
fig.savefig("monq_genotype" + str(nonneut_3) + "_nntrburnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

mod.plot_genotype(locus=nonneut_4, individs=random_individs, edge_color='face', size=0.1, alpha=1) #non-neutral locus from non-neutral data csv file
fig = plt.gcf()
fig.savefig("monq_genotype" + str(nonneut_4) + "_nntrburnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

plt.rcParams["figure.figsize"] = (6,6) 
mod.plot_fitness(individs=random_individs, edge_color='face', size=0.5, alpha=0.5, fit_cmap='seismic')
fig = plt.gcf()
plt.subplots_adjust(right=0.8) 
fig.savefig("monq_fitness_nntrburnin.png", dpi=300, bbox_inches='tight')
plt.close(fig)
  
# After non-neutral runs, check which loci were selected as non-neutral ones and modify below plot codes accordingly.
 
# This is the real, main runs
count = 60
while count < 150:
  mod.walk(10, 'main')
  
  plt.rcParams["figure.figsize"] = (9,4)
  
  random_individs = mod.get_random_individs(n=10000, spp=0)
  
  mod.plot(lyr=0, spp=0, individs=random_individs, color='none', edge_color='black', size=0.1, alpha=0.5)
  fig = plt.gcf()
  fig.savefig("monq_mod_" + str(count) + ".png", dpi=300, bbox_inches='tight')
  plt.close(fig)

  mod.plot_phenotype(individs=random_individs, edge_color='face', size=0.1, alpha=1)
  fig = plt.gcf()
  fig.savefig("monq_phenotype_" + str(count) + ".png", dpi=300, bbox_inches='tight')
  plt.close(fig)

  #density raster is not being drawn
  #mod.plot_density(spp=0, individs=random_individs, size=0.1)
  #fig = plt.gcf()
  #fig.savefig("monq_density_" + str(count) + ".png", dpi=300, bbox_inches='tight')
  #plt.close(fig)

  mod.plot_genetic_PCA(individs=random_individs, edge_color='face', size=0.1, alpha=1)
  fig = plt.gcf()
  fig.savefig("monq_pca_" + str(count) + ".png", dpi=300, bbox_inches='tight')
  plt.close(fig)

  mod.plot_allele_frequencies(spp=0, color='red')
  fig = plt.gcf()
  fig.savefig("monq_af_" + str(count) + ".png", dpi=300, bbox_inches='tight')
  plt.close(fig)

  mod.plot_hist_fitness()
  fig = plt.gcf()
  fig.savefig("monq_hist_ft_" + str(count) + ".png", dpi=300, bbox_inches='tight')
  plt.close(fig)

  mod.plot_movement_surface(style='hist')
  fig = plt.gcf()
  fig.savefig("monq_movement_sf_" + str(count) + ".png", dpi=300, bbox_inches='tight')
  plt.close(fig)

  mod.plot_dispersal_surface(style='hist')
  fig = plt.gcf()
  fig.savefig("monq_dispersal_sf_" + str(count) + ".png", dpi=300, bbox_inches='tight')
  plt.close(fig)

  mod.plot_dispersal(n_individs=10000, include_start_points=False, size=0.1)
  fig = plt.gcf()
  fig.savefig("monq_dispersal_" + str(count) + ".png", dpi=300, bbox_inches='tight')
  plt.close(fig)

  mod.plot_movement(n_individs=10000, include_start_points=False, size=0.1)
  fig = plt.gcf()
  fig.savefig("monq_movement_" + str(count) + ".png", dpi=300, bbox_inches='tight')
  plt.close(fig)

  mod.plot_demographic_pyramid()
  fig = plt.gcf()
  fig.savefig("monq_pyramid_" + str(count) + ".png", dpi=300, bbox_inches='tight')
  plt.close(fig)

  mod.plot_pop_growth(expected=True, actual=True, expected_color='red', actual_color='blue')
  fig = plt.gcf()
  fig.savefig("monq_growth_" + str(count) + ".png", dpi=300, bbox_inches='tight')
  plt.close(fig)

  mod.plot_genotype(locus=nonneut_1, individs=random_individs, edge_color='face', size=0.1, alpha=1) #non-neutral locus from non-neutral data csv file
  fig = plt.gcf()
  fig.savefig("monq_genotype" + str(nonneut_1) + "_" + str(count) + ".png", dpi=300, bbox_inches='tight')
  plt.close(fig)

  mod.plot_genotype(locus=nonneut_2, individs=random_individs, edge_color='face', size=0.1, alpha=1) #non-neutral locus from non-neutral data csv file
  fig = plt.gcf()
  fig.savefig("monq_genotype" + str(nonneut_2) + "_" + str(count) + ".png", dpi=300, bbox_inches='tight')
  plt.close(fig)

  mod.plot_genotype(locus=nonneut_3, individs=random_individs, edge_color='face', size=0.1, alpha=1) #non-neutral locus from non-neutral data csv file
  fig = plt.gcf()
  fig.savefig("monq_genotype" + str(nonneut_3) + "_" + str(count) + ".png", dpi=300, bbox_inches='tight')
  plt.close(fig)

  mod.plot_genotype(locus=nonneut_4, individs=random_individs, edge_color='face', size=0.1, alpha=1) #non-neutral locus from non-neutral data csv file
  fig = plt.gcf()
  fig.savefig("monq_genotype" + str(nonneut_4) + "_" + str(count) + ".png", dpi=300, bbox_inches='tight')
  plt.close(fig)

  plt.rcParams["figure.figsize"] = (6,6) 
  mod.plot_fitness(individs=random_individs, edge_color='face', size=0.5, alpha=0.5, fit_cmap='seismic')
  fig = plt.gcf()
  plt.subplots_adjust(right=0.8) 
  fig.savefig("monq_fitness_" + str(count) + ".png", dpi=300, bbox_inches='tight')
  plt.close(fig)
  
  count += 10
