#!/usr/bin/env python2
# This script can be used to compute 1D and 2D PMFs, as well as the free energy and expectation
# values of predefined clusters.
# 1/30/14 Extending to be able to perform free energy perturbation calculations.
# This script uses the pymbar package: https://github.com/choderalab/pymbar
# To learn how to use the script: $ python ./compute-pmf.py -help
import numpy
import math
import sys
import os
# numpy.set_printoptions(threshold=numpy.nan) # output full matrices when printing
numpy.set_printoptions(threshold=sys.maxsize) # output full matrices when printing


import pickle # used to save files so that the initialization calculations don't have to be repeated
#import pymbar # used to do the computations of free energy
#sys.path.append('C:\\Users\\Dell\\Anaconda\\Lib\\site-packages\\pymbar\\')
#sys.path.append('/home/hht1/programs/anaconda/pkgs/pymbar-3.0.0.beta2-np19py27_1/lib/python2.7/site-packages/pymbar-3.0.0.dev0-py2.7-linux-x86_64.egg/pymbar')
#sys.path.append('/home/hht1/programs/anaconda/lib/python2.7/site-packages/pymbar-3.0.0.dev0-py2.7-linux-x86_64.egg/pymbar/')
#sys.path.append('/home/hht1/CODES/pymbar-2.1.0-beta/build/lib.linux-x86_64-2.7/pymbar')
#sys.path.append('/home/wl45/pymbar/build/lib.linux-ppc64-2.7')
# sys.path.append('/home/wl45/python/lib/python2.7/site-packages/pymbar-3.0.0.dev0-py2.7-linux-ppc64.egg/pymbar')
import pymbar
from pymbar import timeseries # used to subsample data so that the samples are uncorrelated
os.system("echo 'Time' > time.info")
# Flags
debug_flag = False
check_for_pickle_files = True

# Command line arguments
if len(sys.argv) == 1:
    print "#####################################################################################"
    print "# If you don't enter any command line arguments, the script will run with defaults. #"
    print "# You can use $ python ./compute-pmf.py -help to see the command line arguments.    #"
    print "#####################################################################################"

# Set default command line arguments
submit_to_cluster = True          # flag for submitting calculation to a cluster (eliminates interactivity)
precision_threshold = 1e-15       # threshold to determine whether or not to attempt to compute expectation values
subsample_trajectories = False      # flag for using pymbar's built-in subsampling features
compute_per_bin_quantities = False  # flag for computing per-bin expectation values
kB = 1.381e-23 * 6.022e23 / 4184.0  # Boltzmann constant in kcal/mol/K
ndim = 2                           # dimension of the pmf, can be 1 or 2
metadata_file = 'metadatafile'         # name of the metadata file
biasing_variable_column = 5        # column in the trajectory file that contains the biasing variable information
energy_column = 4                  # column in the trajectory file that contains the total (unbiased) potential energy
pmf_variable_column_1 = 5          # column in the trajectory file that contains the first pmf variable
pmf_variable_column_2 = 6          # column in the trajectory file that contians the second pmf variable
nbins1 = 40                        # number of bins for the first pmf variable
nbins2 = 30                        # number of bins for the second pmf variable
start_temperature = 250            # temperature at which to start the free energy calculations
end_temperature = 650              # temperature at which to stop the free energy calculations
force = 0.0
temperature_increment = 10         # how often (in degrees) to compute the free energy
N_samples = 2000                   # number of correlated samples per simulation
expectation_columns = []           # an array of column numbers for which to compute expectation values
expectation_files = []             # an array containing the names of the files that the expectation value data should be read from
pmf_variable_column_1_file = 'DEFAULT'
pmf_variable_column_2_file = 'DEFAULT'
expectation_data = []              # an array of numpy arrays to hold data used to compute expectation values
fep_columns = []                   # an array of column numbers for which to compute free energy perturbations
perturb_from_columns = []          # an array of column numbers which tell MBAR how to perturb away from a particular energy column
fep_data = []                      # an array of numpy arrays to hold data used to compute free energy perturbations
cluster_binning = False            # True if you are using pre-assigned clusters to bin the data
cluster_bin_map = []               # Used to map cluster indices
nbiases = 1                        # Number of biases applied during the simulations
nperturbations = 0                 # Number of perturbed Hamiltonians to evaluate
biasing_variable_columns = []      # List of biasing variable columns in data file
biasing_variable_columns.append(biasing_variable_column)
specified_a_biasing_variable_column = False
perturbation_columns = []          # List of columns in data file with perturbed Hamiltonian energies
biasing_variable_kt = []           # List of numpy arrays containing biasing variables for all snapshots
biasing_variable_kn = []           # List of numpy arrays containing biasing variables for all snapshots (used to store subsampled data)
U_kt = []                          # List of numpy arrays containing perturbed energies for all snapshots
U_kn = []                          # List of numpy arrays containing perturbed energies for all snapshots (used to store subsampled data)
u_kn = []                          # List of numpy arrays containing reduced perturbed energies for all snapshots (used to store subsampled data)

# Read command line arguments
for arg in range(len(sys.argv)):
  if arg % 2 == 0:
    continue
  if sys.argv[arg] == "-kb":
    kB = float(sys.argv[arg+1])
  elif sys.argv[arg] == "-m":
    metadata_file = sys.argv[arg+1]
  elif sys.argv[arg] == "-ev":
    ev_argument = sys.argv[arg+1]
    if ':' in ev_argument:
      ev_file, ev_columns = ev_argument.split(':')
    else:
      ev_file = 'DEFAULT'
      ev_columns = ev_argument
    if '-' in ev_columns:
      ev_columns = ev_columns.split('-')
      starting_column = int(ev_columns[0])
      ending_column = int(ev_columns[1])
      for i in range(starting_column, ending_column+1):
        expectation_columns.append(i)
        expectation_files.append(ev_file)
    else:
      expectation_columns.append(int(ev_columns))
      expectation_files.append(ev_file)
  elif sys.argv[arg] == "-fep":
    argument_to_split = sys.argv[arg+1]
    fep_column, perturb_from_column = argument_to_split.split(",")
    fep_columns.append(int(fep_column))
    perturb_from_columns.append(int(perturb_from_column))
  elif sys.argv[arg] == "-b":
    # Clear the default biasing variable column if you explicitly specify a column
    if specified_a_biasing_variable_column == False:
      specified_a_biasing_variable_column = True
      biasing_variable_columns = []
    biasing_variable_column = int(sys.argv[arg+1])
    biasing_variable_columns.append(biasing_variable_column)
  elif sys.argv[arg] == "-e":
    energy_column = int(sys.argv[arg+1])
  elif sys.argv[arg] == "-d":
    ndim = int(sys.argv[arg+1])
    if ndim !=1 and ndim !=2:
      print "ndim must be either 1 or 2"
      sys.exit()
  elif sys.argv[arg] == "-v1":
    pmf_variable_column_1 = sys.argv[arg+1]
    if ':' in pmf_variable_column_1:
      pmf_variable_column_1_file, pmf_variable_column_1 = pmf_variable_column_1.split(':')
      pmf_variable_column_1 = int(pmf_variable_column_1)
    else:
      pmf_variable_column_1_file = 'DEFAULT'
      pmf_variable_column_1 = int(pmf_variable_column_1)
  elif sys.argv[arg] == "-v1n":
    nbins1 = int(sys.argv[arg+1])
  elif sys.argv[arg] == "-v2":
    pmf_variable_column_2 = sys.argv[arg+1]
    if ':' in pmf_variable_column_2:
      pmf_variable_column_2_file, pmf_variable_column_2 = pmf_variable_column_2.split(':')
      pmf_variable_column_2 = int(pmf_variable_column_2)
    else:
      pmf_variable_column_2_file = 'DEFAULT'
      pmf_variable_column_2 = int(pmf_variable_column_2)
  elif sys.argv[arg] == "-v2n":
    nbins2 = int(sys.argv[arg+1])
  elif sys.argv[arg] == "-st":
    start_temperature = int(sys.argv[arg+1])
  elif sys.argv[arg] == "-et":
    end_temperature = int(sys.argv[arg+1])
  elif sys.argv[arg] == "-f":
    force = float(sys.argv[arg+1])
  elif sys.argv[arg] == "-ti":
    temperature_increment = int(sys.argv[arg+1])
  elif sys.argv[arg] == "-nsamples":
    N_samples = int(sys.argv[arg+1])
  elif sys.argv[arg] == "-ss":
    ss = sys.argv[arg+1]
    if ss != 'y' and ss != 'n':
      print "subsample_trajectories must be either y or n"
      sys.exit()
    if ss == 'y':
      subsample_trajectories = True
    elif ss == 'n':
      subsample_trajectories = False
  elif sys.argv[arg] == "-pb":
    pb = sys.argv[arg+1]
    if pb != 'y' and pb != 'n':
      print "compute_per_bin_quantities must be either y or n"
      sys.exit()
    if pb == 'y':
      compute_per_bin_quantities = True
    elif pb == 'n':
      compute_per_bin_quantities = False
  elif sys.argv[arg] == "-submit":
    submit = sys.argv[arg+1]
    if submit != 'y' and submit != 'n':
      print "submit must be either y or n"
      sys.exit()
    if submit == 'y':
      submit_to_cluster = True
    elif submit == 'n':
      submit_to_cluster = False
  elif sys.argv[arg] == "-clustering":
      if sys.argv[arg+1] == 'y':
          cluster_binning = True
  elif sys.argv[arg] == "-p":
      perturbation_columns.append(int(sys.argv[arg+1]))
  elif sys.argv[arg] == "-help":
    print "Usage: $ python ./compute-pmf.py \\ \n -kb boltzmann_constant (default=%f)\\ \n -m metadatafile (default=%s)\\ \n -b biasing_variable_column (default=%d)\\ \n -e energy_column (default=%d)\\ \n -d pmf_dimension (default=%d)\\ \n -v1 pmf_variable_column_1 (default=%d)\\ \n -v1n pmf_variable_nbins_1 (default=%d)\\ \n -v2 pmf_variable_column_2 (default=%d)\\ \n -v2n pmf_variable_nbins_2 (default=%d)\\ \n -st start_temperature (default=%d)\\ \n -et end_temperature (default=%d)\\ \n -f force (default=%f)\\ \n -ti temperature_increment (default=%d)\\ \n -nsamples N_samples (default=%d)\\ \n -ev column_number_to_compute_expection_value (no default) \\ \n -pb compute_per_bin_quantities (default=%s)\\ \n -ss subsample_trajectories (default=%s)\\ \n -submit submit_to_cluster (default=%s)\\ \n -clustering cluster_binning (default=%s)\\ \n -p perturbation_column (default=None)\\ \n" % (kB, metadata_file, biasing_variable_column, energy_column, ndim, pmf_variable_column_1, nbins1, pmf_variable_column_2, nbins2, start_temperature, end_temperature, force, temperature_increment, N_samples, compute_per_bin_quantities, subsample_trajectories, submit_to_cluster, cluster_binning)
    print "All command line arguments are optional and can be specified in an arbitrary order on a single line. To change the defaults, use a text editor to edit the script."
    sys.exit()
  else:
    print "I don't recognize %s as an argument. Use -help to see the proper usage." % sys.argv[arg]
    sys.exit()

# exec(open("config.py").read())

nbiases = len(biasing_variable_columns)
nperturbations = len(perturbation_columns)

print "Using Boltzmann constant value kB=%f" % kB
print "Using metadata file %s" % metadata_file
print "Getting %d biasing variable(s) data from column(s) %s" % (nbiases, biasing_variable_columns)
print "Performing %d perturbation calculations using data from column(s) %s" % (nperturbations, perturbation_columns)
print "Getting energy data from column %d" % (energy_column)
print "Getting pmf variable 1 data from column %d" % (pmf_variable_column_1)
print "Using %d bins to bin pmf variable 1" % nbins1
if ndim == 2:
  print "Getting pmf variable 2 data from column %d" % (pmf_variable_column_2)
  print "Using %d bins to bin pmf variable 2" % nbins2
print "Starting from a temperature %d K and ending at a temperature %d K and using a temperature increment of %d K." % (start_temperature, end_temperature, temperature_increment)
print "#############################################################"
print "# Warning: you must have exactly %d samples per simulation #" % N_samples
print "#############################################################"
if subsample_trajectories == True:
  print "Subsampling of trajectories is ON."
elif subsample_trajectories == False:
  print "Subsampling of trajectories is OFF."
if len(expectation_columns) > 0:
  print "Calculating expectation values for the following columns: %s \n From files: %s" % (str(expectation_columns), str(expectation_files))
if len(fep_columns) > 0:
  print "Calculating free energy perturbation values for the following columns: %s" % str(fep_columns)
if compute_per_bin_quantities == True:
  print "Computing per bin quantities is ON."
elif compute_per_bin_quantities == False:
  print "Computing per bin quantities is OFF."

if submit_to_cluster == False:
  while 1:
    answer = raw_input("Would you like to continue? (y/n)\n")
    if answer == 'y':
      break
    elif answer == 'n':
      print "Okay. Exiting."
      sys.exit()

  # Check to see if you've already generated an mbar.pkl file and offer to reuse it to save time
  load_pickle = False
  if check_for_pickle_files:
      if os.path.exists('./mbar.pkl'):
          while 1:
              answer = raw_input("mbar.pkl exists. Would you like to load the file (and skip the initialization of mbar)? y/n\n")
              if answer == 'y':
                  load_pickle = True
                  break
              elif answer == 'n':
                  load_pickle = False
                  break
elif submit_to_cluster == True:
  # Check to see if you've already generated an mbar.pkl file and reuse it to save time
  load_pickle = False
  if check_for_pickle_files:
      if os.path.exists('./mbar.pkl'):
          load_pickle = True
          print "Using the existing mbar.pkl; skipping the initialization of the calculation."
      else:
          load_pickle = False
          print "Could not find an existing mbar.pkl; starting the initialization of the calculation."

# If not loading pickle files, start the intitialization calculation
if load_pickle == False:
    # Lists
    files = []              # list of trajectory files
    temperatures = []       # list of sampling temperatures
    biasing_strengths = []  # list of biasing strengths
    biasing_values = []     # list of biasing center values
    for i in range(nbiases):
        biasing_strengths.append([])
        biasing_values.append([])

    # Read metadata file
    metadata = open(metadata_file, 'r')
    for line in metadata:
        column_index = 0
        if line[0] == '#':
            continue
        line = line.split()
        files.append(line[column_index])
        column_index += 1
        temperatures.append(float(line[column_index]))
        column_index += 1
        for i in range(nbiases):
            biasing_strengths[i].append(float(line[column_index]))
            column_index += 1
            biasing_values[i].append(float(line[column_index]))
            column_index += 1
    metadata.close()

    if debug_flag:
        print "Information read from metadata file:"
        print "Files:"
        print files
        print "Temperatures:"
        print temperatures
        print "Biasing strengths:"
        print biasing_strengths
        print "Biasing values:"
        print biasing_values

    # Parameters
    K = len(files)  # number of simulations

    # Allocate storage for simulation data
    N_k = numpy.zeros([K], numpy.int32)  # N_k[k] is the number of snapshots from umbrella simulation k
    for i in range(nbiases):
        biasing_variable_kt.append(numpy.zeros([K,N_samples], numpy.float64))
    for i in range(nperturbations+1):
        U_kt.append(numpy.zeros([K,N_samples], numpy.float64))
    pmf_variable_kt_1 = numpy.zeros([K,N_samples], numpy.float64)  # pmf_variable_kt_1[k,t] is the value of the first pmf variable from snapshot t of simulation k
    if ndim == 2:
        pmf_variable_kt_2 = numpy.zeros([K,N_samples], numpy.float64)  # pmf_variable_kt_2[k,t] is the value of the second pmf variable from snapshot t of simulation k
    cluster_bin_kt = -1*numpy.ones([K,N_samples], numpy.int32)  # cluster_bin_kt[k,t] is the bin of snapshot t of umbrella simulation k
    if (len(expectation_columns) > 0):
        for i in range(len(expectation_columns)):
            expectation_data.append(numpy.zeros([K,N_samples], numpy.float64))
    if (len(fep_columns) > 0):
        for i in range(len(fep_columns)):
            fep_data.append(numpy.zeros([K,N_samples], numpy.float64))

    # Read the correlated simulation data.
    print "Reading simulation data from default files..."
    k =0 # simuluation index
    for trajectory_file in files:
        print "Reading file %s ..." % trajectory_file
        kT = kB * temperatures[k]   # thermal energy
        beta = 1.0 / kT             # inverse temperature
        t = 0                       # time step index
        for line in open(trajectory_file, 'r'):
            if line[0] == '#':
                continue
            line = line.split()
            # Store order parameters.
            for i in range(nbiases):
                # print i
                # print biasing_variable_columns[i]
                # print biasing_variable_columns
                # print float(line[biasing_variable_columns[i]-1])
                # print k, t
                biasing_variable_kt[i][k,t] = float(line[biasing_variable_columns[i]-1])
            if pmf_variable_column_1_file == 'DEFAULT':
              pmf_variable_kt_1[k,t] = float(line[pmf_variable_column_1-1])
            if ndim == 2:
              if pmf_variable_column_2_file == 'DEFAULT':
                pmf_variable_kt_2[k,t] = float(line[pmf_variable_column_2-1])
            U_kt[0][k,t] = float(line[energy_column-1])
            for i in range(nperturbations):
                U_kt[i+1][k,t] = float(line[perturbation_columns[i]-1])
            if len(expectation_columns) > 0:
                for i in range(len(expectation_columns)):
                    if expectation_files[i] == 'DEFAULT':
                      expectation_data[i][k,t] = float(line[expectation_columns[i]-1])
            if len(fep_columns) > 0:
                for i in range(len(fep_columns)):
                    fep_data[i][k,t] = float(line[perturb_from_columns[i]-1])-float(line[fep_columns[i]-1])
            t += 1 # increment time step index
        k += 1     # increment simulation index

    # Read the correlated simulation data from other files
    print "Reading simulation data from other files..."
    k = 0 # simuluation index
    expectation_files_set = list(set(expectation_files))
    for trajectory_file in files:
        directory = os.path.dirname(os.path.realpath(trajectory_file))
        print "Reading extra files for the trajectory file %s ..." % trajectory_file
        if pmf_variable_column_1_file != 'DEFAULT':
          t = 0                       # time step index
          for line_index, line in enumerate(open(os.path.join(directory,pmf_variable_column_1_file), 'r')):
              if line_index == 0:
                  print "WARNING: Skipping first line of %s by default." % pmf_variable_column_1_file
                  continue
              line = line.split()
              # Store order parameters.
              pmf_variable_kt_1[k,t] = float(line[pmf_variable_column_1-1])
              t += 1 # increment time step index
        if pmf_variable_column_2_file != 'DEFAULT':
          t = 0                       # time step index
          for line_index, line in enumerate(open(os.path.join(directory,pmf_variable_column_2_file), 'r')):
              if line_index == 0:
                  print "WARNING: Skipping first line of %s by default." % pmf_variable_column_2_file
                  continue
              line = line.split()
              # Store order parameters.
              pmf_variable_kt_2[k,t] = float(line[pmf_variable_column_2-1])
              t += 1 # increment time step index
        for expectation_file in expectation_files_set:
          t = 0                       # time step index
          if expectation_file == "DEFAULT":
              break
          for line_index, line in enumerate(open(os.path.join(directory,expectation_file), 'r')):
              if line_index == 0:
                  print "WARNING: Skipping first line of %s by default." % expectation_file
                  continue
              line = line.split()
              # Store order parameters.
              if len(expectation_columns) > 0:
                  for i in range(len(expectation_columns)):
                      if expectation_file == expectation_files[i]:
                        expectation_data[i][k,t] = float(line[expectation_columns[i]-1])
              t += 1 # increment time step index
        k += 1     # increment simulation index

    # If using clusters, read the correlated cluster bin indices
    if cluster_binning:
    # k == trajectory index
        k = 0
        for trajectory_file in files:
            cluster_file = "/".join(trajectory_file.split("/")[:-1])+"/cluster.dat"
            if not os.path.exists(cluster_file):
                print "%s does not exist. Exiting." % cluster_file
                sys.exit()
            f = open(cluster_file, "r")
            # n == snapshot index
            n = 0
            for line in f:
                line = line.split()
                if line[0] == "#": continue
                cluster_bin_kt[k,n] = int(line[0])
                n += 1
            f.close()
            k += 1

    # Subsample data.
    print "Subsampling data..."
    for i in range(nbiases):
        biasing_variable_kn.append(numpy.zeros([K,N_samples], numpy.float64)) # biasing_variable_kn[k,n] is the value of the biasing variable from snapshot n of simulation k
    for i in range(nperturbations+1):
        U_kn.append(numpy.zeros([K,N_samples], numpy.float64)) # biasing_variable_kn[k,n] is the value of the biasing variable from snapshot n of simulation k
    if not cluster_binning:
        pmf_variable_kn_1 = numpy.zeros([K,N_samples], numpy.float64) # pmf_variable_kn_1[k,n] is the value of the first pmf variable from snapshot n of simulation k
        if ndim == 2:
            pmf_variable_kn_2 = numpy.zeros([K,N_samples], numpy.float64) # pmf_variable_kn_2[k,n] is the value of the second pmf variable from snapshot n of simulation k
    else:
        cluster_bin_kn = -1*numpy.ones([K,N_samples], numpy.int32) # cluster_bin_kn[k,n] is the cluster bin index of snapshot n of umbrella simulation k
    N_k = numpy.zeros([K], numpy.int32) # N_k[k] is the number of uncorrelated samples from simulation index k
    reduced_expectation_data = []
    if len(expectation_columns) > 0:
        for i in range(len(expectation_columns)):
            reduced_expectation_data.append(numpy.zeros([K,N_samples], numpy.float64))
    reduced_fep_data = []
    if len(fep_columns) > 0:
        for i in range(len(fep_columns)):
            reduced_fep_data.append(numpy.zeros([K,N_samples], numpy.float64))
    for k in range(K):
        # Extract timeseries.
        A_t = biasing_variable_kt[0][k,:]
        # Compute statistical inefficiency.
        try:
            g = timeseries.statisticalInefficiency(A_t)
        except Exception as e:
            print str(e)
            print A_t

        # Subsample data.
        if subsample_trajectories:
            indices = timeseries.subsampleCorrelatedData(A_t, g=g)
        else:
            indices = timeseries.subsampleCorrelatedData(A_t, g=1)
        N = len(indices) # number of uncorrelated samples
        print "k = %5d : g = %.1f, N = %d" % (k, g, N)
        for i in range(nbiases):
            biasing_variable_kn[i][k,0:N] = biasing_variable_kt[i][k,indices]
        for i in range(nperturbations+1):
            U_kn[i][k,0:N] = U_kt[i][k,indices]
        if not cluster_binning:
            pmf_variable_kn_1[k,0:N] = pmf_variable_kt_1[k,indices]
            if ndim == 2:
                pmf_variable_kn_2[k,0:N] = pmf_variable_kt_2[k,indices]
        if cluster_binning:
            cluster_bin_kn[k,0:N] = cluster_bin_kt[k,indices]
        if len(expectation_columns) > 0:
            for i in range(len(expectation_columns)):
                reduced_expectation_data[i][k,0:N] = expectation_data[i][k,indices]
        if len(fep_columns) > 0:
            for i in range(len(fep_columns)):
                reduced_fep_data[i][k,0:N] = fep_data[i][k,indices]
        N_k[k] = N # number of uncorrelated samples

    # Find unique cluster indices
    if cluster_binning:
        cluster_bin_map = numpy.unique(cluster_bin_kn)
        cluster_bin_map = cluster_bin_map[cluster_bin_map > 0]
        print "Clusters found: " + str(cluster_bin_map)

    # Compute reduced potentials from all simulations in all thermodynamic states.
    print "Computing reduced potentials..."
    u_kln = numpy.zeros([K,K,N_samples], numpy.float64) # u_kln[k,l,n] is reduced biased potential of uncorrelated snapshot n from simulation k in thermodynamic state l
    for k in range(K):
        N = N_k[k] # number of uncorrelated snapshots
        for l in range(K):
            # Compute reduced potential in all other temperatures and biasing potentials indexed by l.
            kT = kB * temperatures[l] # thermal energy
            beta = 1.0 / kT           # inverse temperature
            u_kln[k,l,0:N] = beta * (U_kn[0][k,0:N])
            for i in range(nbiases):
                U_bias = (biasing_strengths[i][l]/2.0) * (biasing_variable_kn[i][k,0:N] - biasing_values[i][l])**2 + (biasing_variable_kn[i][k,0:N]-25.1)*force   # biasing potential for this sample
                # U_bias = (biasing_variable_kn[i][l,0:N]-25.1)*force   # biasing potential for this sample
                # print "------"
                # print biasing_strengths[i][l]/2.0
                # print biasing_values[i][l]
                # print biasing_variable_kn[i][k,0:10]
                # print N
                # print biasing_variable_kn[i][k,0:10] - biasing_values[i][l]
                # print (biasing_strengths[i][l]/2.0) * (biasing_variable_kn[i][k,0:10] - biasing_values[i][l])**2
                # print "------"
                # print (biasing_variable_kn[i][k,0:10]-25.1)*force
                # print biasing_variable_kn[i][k,0:10]
                # print kT
                # print "------"

                u_kln[k,l,0:N] += beta * (U_bias)
            # exit()
    # Bin data
    print "Binning data..."
    bin_kn = -1 * numpy.ones([K,N_samples], numpy.int32) # bin_kn[k,n] is bin index of sample n from simulation k; otherwise -1
    if not cluster_binning:
        # Construct PMF bins in pmf_variable
        pmf_variable_min_1 = pmf_variable_kt_1.min()
        pmf_variable_max_1 = pmf_variable_kt_1.max()
        pmf_variable_bin_width_1 = (pmf_variable_max_1 - pmf_variable_min_1) / nbins1
        if ndim == 2:
            pmf_variable_min_2 = pmf_variable_kt_2.min()
            pmf_variable_max_2 = pmf_variable_kt_2.max()
            pmf_variable_bin_width_2 = (pmf_variable_max_2 - pmf_variable_min_2) / nbins2

        # Bin data
        bin_centers = list() # bin_centers[i] is the center of bin i
        if ndim == 1:
            for i in range(nbins1):
                bin_center_1 = pmf_variable_min_1 + pmf_variable_bin_width_1 * (i + 0.5)
                bin_centers.append( bin_center_1 )
        elif ndim == 2:
            for i in range(nbins1):
                for j in range(nbins2):
                    bin_center_1 = pmf_variable_min_1 + pmf_variable_bin_width_1 * (i + 0.5)
                    bin_center_2 = pmf_variable_min_2 + pmf_variable_bin_width_2 * (j + 0.5)
                    bin_centers.append( (bin_center_1, bin_center_2) )
        # print("bin_center--------------")
        # print(bin_centers)
        if ndim == 1:
            for k in range(K):
                N = N_k[k]
                # Compute bin assignment for 1D case.
                bin_kn[k,0:N] = numpy.int32((pmf_variable_kn_1[k,0:N] - pmf_variable_min_1) / pmf_variable_bin_width_1) # pmf bin index 0..(nbins-1)
        elif ndim == 2:
            for k in range(K):
                N = N_k[k]
                # Compute bin assignment for 2D case.
                bin_1 = numpy.int32((pmf_variable_kn_1[k,0:N] - pmf_variable_min_1) / pmf_variable_bin_width_1) # pmf bin index 0..(nbins-1)
                bin_2 = numpy.int32((pmf_variable_kn_2[k,0:N] - pmf_variable_min_2) / pmf_variable_bin_width_2) # pmf bin index 0..(nbins-1)
                bin_kn[k,0:N] = bin_1 * nbins2 + bin_2
    else:
        for k in range(K):
            N = N_k[k]
            for n in range(N):
                bin_kn[k,n] = cluster_bin_map[numpy.where(cluster_bin_map == cluster_bin_kn[k,n])[0][0]]

    # Make sure all bins are populated.
    print "Counting number of samples per bin..."
    if cluster_binning:
        nbins = cluster_bin_map.size
        bin_counts = numpy.zeros([nbins], numpy.int32)
        for i in range(nbins):
            bin_counts[i] = (bin_kn == cluster_bin_map[i]).sum()
            # I set up a cutoff to simplify the view
    else:
        if ndim == 1:
            nbins = nbins1
            bin_counts = numpy.zeros([nbins], numpy.int32)
        if ndim == 2:
            nbins = nbins1*nbins2
            bin_counts = numpy.zeros([nbins], numpy.int32)
        for i in range(nbins):
            bin_counts[i] = (bin_kn == i).sum()
    print "Number of samples in each bin:"
    print bin_counts

    # Map to reduced bin counts to eliminate empty bins.
    nonempty_bins = numpy.where(bin_counts > 0)[0]
    # nonempty_bins = numpy.where(bin_counts > 10)[0]
    # nonempty_bins = numpy.where(bin_counts > 2)[0]
    # nonempty_bins = numpy.where(bin_counts > 5)[0]
    print "Non-empty bins:"
    print nonempty_bins
    nreduced_bins = len(nonempty_bins)
    reduced_bin_kn = -1 * numpy.ones([K,N_samples], numpy.int32) # reduced_bin_kn[k,n] is bin index of sample n from simulation k; otherwise -1
    bin_map = [] # save original bin number so that you can map the free energy back to its original place
    if cluster_binning:
        for (i,j) in enumerate(nonempty_bins):
            bin_map.append(j)
            indices = numpy.where(bin_kn == cluster_bin_map[j])
            reduced_bin_kn[indices] = i
    else:
        for (i,j) in enumerate(nonempty_bins):
            bin_map.append(j)
            indices = numpy.where(bin_kn == j)
            reduced_bin_kn[indices] = i

    # Build in_this_bin arrays, used for computing free energies and expectation values
    in_this_bin = []
    for i in range (nreduced_bins):
        in_this_bin.append(numpy.where(reduced_bin_kn == i, 1.0, 0.0))

    # Initialize MBAR.
    print "Running MBAR..."
    mbar = pymbar.MBAR(u_kln, N_k, verbose=True)
    # Pickle files for later use
    pickle.dump(mbar, open("mbar.pkl", "wb"))
    pickle.dump(U_kn, open("Ukn.pkl", "wb"))
    pickle.dump(nperturbations, open("nperturbations.pkl", "wb"))
    pickle.dump(reduced_bin_kn, open("reducedbinkn.pkl", "wb"))
    pickle.dump(nreduced_bins, open("nreducedbins.pkl", "wb"))
    pickle.dump(K, open("k.pkl", "wb"))
    pickle.dump(N_k, open("nk.pkl", "wb"))
    pickle.dump(bin_map, open("binmap.pkl", "wb"))
    pickle.dump(temperatures, open("temperatures.pkl", "wb"))
    pickle.dump(bin_counts, open("bincounts.pkl", "wb"))
    pickle.dump(reduced_expectation_data, open("reducedexpectationdata.pkl", "wb"))
    pickle.dump(reduced_fep_data, open("reducedfepdata.pkl", "wb"))
    pickle.dump(in_this_bin, open("inthisbin.pkl", "wb"))
    pickle.dump(fep_columns, open("fepcolumns.pkl", "wb"))
    pickle.dump(perturb_from_columns, open("perturbfromcolumns.pkl", "wb"))
    if not cluster_binning:
        pickle.dump(bin_centers, open("bincenters.pkl", "wb"))
        pickle.dump(pmf_variable_kn_1, open("pmfvariablekn1.pkl", "wb"))
        if ndim == 2:
          pickle.dump(pmf_variable_kn_2, open("pmfvariablekn2.pkl", "wb"))
    else:
        pickle.dump(cluster_bin_map, open("clusterbinmap.pkl", "wb"))

else:
    # Load files to save time
    mbar = pickle.load(open("mbar.pkl", "rb"))
    U_kn = pickle.load(open("Ukn.pkl", "rb"))
    nperturbations = pickle.load(open("nperturbations.pkl", "rb"))
    reduced_bin_kn = pickle.load(open("reducedbinkn.pkl", "rb"))
    nreduced_bins = pickle.load(open("nreducedbins.pkl", "rb"))
    K = pickle.load(open("k.pkl", "rb"))
    N_k = pickle.load(open("nk.pkl", "rb"))
    bin_map = pickle.load(open("binmap.pkl", "rb"))
    temperatures = pickle.load(open("temperatures.pkl", "rb"))
    bin_counts = pickle.load(open("bincounts.pkl", "rb"))
    reduced_expectation_data = pickle.load(open("reducedexpectationdata.pkl", "rb"))
    reduced_fep_data = pickle.load(open("reducedfepdata.pkl", "rb"))
    in_this_bin = pickle.load(open("inthisbin.pkl", "rb"))
    fep_columns = pickle.load(open("fepcolumns.pkl", "rb"))
    perturb_from_columns = pickle.load(open("perturbfromcolumns.pkl", "rb"))
    if not cluster_binning:
        bin_centers = pickle.load(open("bincenters.pkl", "rb"))
        pmf_variable_kn_1 = pickle.load(open("pmfvariablekn1.pkl", "rb"))
        if ndim == 2:
          pmf_variable_kn_2 = pickle.load(open("pmfvariablekn2.pkl", "rb"))
    else:
        cluster_bin_map = pickle.load(open("clusterbinmap.pkl", "rb"))

# Loop over original energy and any perturbed energies, if available
for perturbation_index in range(nperturbations+1):
    # If this is the calculation using the original energies, do not attach a prefix to the files
    if perturbation_index == 0:
        file_prefix = ""
    # If this is a calculation using perturbed energies, attach a prefix
    else:
        file_prefix = "perturbation-%s-" % str(perturbation_index)
        print " ############################################################# "
        print " # Computing pmfs and expectation values for perturbation %s #" % str(perturbation_index)
        print " ############################################################# "

    cv_file = file_prefix + "cv-%s-%s-%s.dat" % (str(start_temperature), str(end_temperature), str(temperature_increment))
    cv_output_file = open(cv_file, 'w')
    cv_output_file.write("# %10s %10s\n" % ('temperature', 'cv'))

    expectation_file = file_prefix + "ev-%s-%s-%s.dat" % (str(start_temperature), str(end_temperature), str(temperature_increment))
    ev_output_file = open(expectation_file, 'w')
    ev_output_file.write("# %10s %10s\n" % ('temperature', 'expectation values'))

    fep_file = file_prefix + "fep-%s-%s-%s.dat" % (str(start_temperature), str(end_temperature), str(temperature_increment))
    fep_output_file = open(fep_file, 'w')
    fep_output_file.write("# %10s %10s\n" % ('temperature', 'fep values'))

    # Compute perturbed reduced potential at each temperature of interest in the absence of a biasing potential.
    for target_temperature in range(start_temperature, end_temperature+1, temperature_increment):
        pmf_file = file_prefix + "pmf-" + str(target_temperature) + ".dat"
        pmf_output_file = open(pmf_file, 'w')
        print "Computing perturbed reduced potential at %d K ..." % target_temperature
        u_kn.append(numpy.zeros([K,N_samples], numpy.float64)) # u_kn[k,n] is the unbiased reduced potential energy of snapshot n of umbrella simulation k at conditions of interest
        kT = kB * target_temperature
        beta = 1.0 / kT # reduced temperature
        for k in range(K):
            N = N_k[k]
            u_kn[perturbation_index][k,0:N] = beta * U_kn[perturbation_index][k,0:N] # unbiased reduced potential at desired temperature

        print "Computing PMF..."
        # Compute PMF in unbiased potential (in units of kT).
        f_i = numpy.zeros(nreduced_bins, numpy.float64)
        df_i = numpy.zeros(nreduced_bins, numpy.float64)
        bin_expectation = numpy.zeros(nreduced_bins, numpy.float64)
        for i in range (nreduced_bins):
            bin_expectation[i] = mbar.computePerturbedExpectation(u_kn[perturbation_index], in_this_bin[i], compute_uncertainty=False)[0]
            if bin_expectation[i] > precision_threshold:
                f_i[i] = -numpy.log(bin_expectation[i])
                # print bin_expectation[i]
            else:
                f_i[i] = numpy.nan
                print bin_expectation[i], i
            df_i[i] = 0.0
        for i in range(nreduced_bins):
            f_i[i] -= numpy.nanmin(f_i)

        print "Computing CV..."
        # Compute <E> and <E^2> to get Cv
        E = mbar.computePerturbedExpectation(u_kn[perturbation_index], U_kn[perturbation_index], compute_uncertainty=False)[0]
        Esquared = mbar.computePerturbedExpectation(u_kn[perturbation_index], U_kn[perturbation_index]**2, compute_uncertainty=False)[0]
        heat_capacity = (Esquared-math.pow(E,2))/(kT*target_temperature)
        if debug_flag:
            print E, Esquared, heat_capacity
        cv_output_file.write("%10d %10.3f\n" % (target_temperature, heat_capacity))

        if len(reduced_expectation_data) > 0:
            print "Computing expectation values as a function of temperature..."
            expectation_values = []
            ev_output_file.write("%10d " % target_temperature)
            for i in range(len(reduced_expectation_data)):
                expectation_values.append(mbar.computePerturbedExpectation(u_kn[perturbation_index], reduced_expectation_data[i], compute_uncertainty=False)[0])
                ev_output_file.write("%10.3f" % expectation_values[i])
            ev_output_file.write("\n")

        # if len(reduced_fep_data) > 0:
        #     print "Computing free energy perturbation values as a function of temperature..."
        #     fep_values = []
        #     fep_output_file.write("%10d " % target_temperature)
        #     for i in range(len(reduced_expectation_data)):
        #         fep_values.append(numpy.exp(beta)*mbar.computePerturbedExpectation(u_kn[perturbation_index], reduced_fep_data[i], compute_uncertainty=False)[0])
        #         fep_output_file.write("%10.3f" % fep_values[i])
        #     fep_output_file.write("\n")

        # Compute expectation value of the energy and entropy for each bin
        E_i = numpy.zeros(nreduced_bins, numpy.float64)
        S_i = numpy.zeros(nreduced_bins, numpy.float64)
        if compute_per_bin_quantities:
            print "Computing per-bin expectation values of energy and entropy..."
            for i in range (nreduced_bins):
                # Check to make sure you will actually be able to compute an expectation value
                if bin_expectation[i] > precision_threshold:
                    E = mbar.computePerturbedExpectation(u_kn[perturbation_index], U_kn[perturbation_index]*in_this_bin[i], compute_uncertainty=False)[0]/bin_expectation[i]
                else:
                    E = numpy.nan
                E_i[i] = E/(kB*target_temperature)
                S_i[i] = E_i[i]-f_i[i]

            # Compute other per-bin expectation values
            if len(reduced_expectation_data) > 0:
                print "Computing other per-bin expectation values..."
                ev_pb_file = file_prefix + "evpb-" + str(target_temperature) + ".dat"
                ev_pb_output_file = open(ev_pb_file, 'w')
                ev_pb_output_file.write("# expectation values for reduced (non-empty) bins\n")
                for i in range(nreduced_bins):
                    if cluster_binning:
                        bin_center_1 = cluster_bin_map[i]
                    else:
                        if ndim == 1:
                            bin_center_1 = bin_centers[bin_map[i]]
                        if ndim == 2:
                            bin_center_1 = bin_centers[bin_map[i]][0]
                            bin_center_2 = bin_centers[bin_map[i]][1]
                    if cluster_binning:
                        ev_pb_output_file.write("  %12d %12d " % (bin_map[i], bin_center_1))
                    else:
                        if ndim == 1:
                            ev_pb_output_file.write("  %12d %12.3f " % (bin_map[i], bin_center_1))
                        if ndim == 2:
                            ev_pb_output_file.write("  %12d %12.3f %12.3f " % (bin_map[i], bin_center_1, bin_center_2))
                    # if bin clustering, include pmf information in evpb file
                    if cluster_binning:
                        ev_pb_output_file.write("%12.3f %12.3f %12.3f %12.3f " % (f_i[i], df_i[i], E_i[i], S_i[i]))
                    # compute each expectation value in list
                    for j in range(len(reduced_expectation_data)):
                        data_in_this_bin = reduced_expectation_data[j]*in_this_bin[i]
                        # Check to make sure you will actually be able to compute an expectation value
                        if bin_expectation[i] > precision_threshold:
                            E = mbar.computePerturbedExpectation(u_kn[perturbation_index], data_in_this_bin, compute_uncertainty=False)[0]/bin_expectation[i]
                        else:
                            E = float('nan')
                        ev_pb_output_file.write("%12.3f " % E)
                    ev_pb_output_file.write("\n")
                ev_pb_output_file.close()

            # Compute other per-bin expectation values
            if len(reduced_fep_data) > 0:
                print "Computing per-bin free energy perturbation values..."
                fep_pb_file = file_prefix + "feppb-" + str(target_temperature) + ".dat"
                fep_pb_output_file = open(fep_pb_file, 'w')
                fep_pb_output_file.write("# free energy perturbation values for reduced (non-empty) bins\n")
                fep_values = numpy.zeros([nreduced_bins, len(reduced_fep_data)], numpy.float64)
                for i in range(nreduced_bins):
                    # compute each expectation value in list
                    for j in range(len(reduced_fep_data)):
                        data_in_this_bin = numpy.exp(reduced_fep_data[j])*in_this_bin[i]
                        # Check to make sure you will actually be able to compute an expectation value
                        if bin_expectation[i] > precision_threshold:
                            E = numpy.log(numpy.exp(beta)*mbar.computePerturbedExpectation(u_kn[perturbation_index], data_in_this_bin, compute_uncertainty=False)[0]/bin_expectation[i])
                            if perturb_from_columns[j] == energy_column:
                                E = f_i[i]-E
                            else:
                                # find the column index where fep_columns == perturb_from_columns[j]
                                column_to_perturb_from = fep_columns.index(perturb_from_columns[j])
                                E = fep_values[i][column_to_perturb_from]-E
                            fep_values[i][j] = E
                        else:
                            E = float('nan')
                            fep_values[i][j] = E
                for j in range(len(reduced_fep_data)):
                    min_value = numpy.nanmin(zip(*fep_values)[j])
                    for i in range(nreduced_bins):
                        fep_values[i][j] -= min_value

                for i in range(nreduced_bins):
                    if cluster_binning:
                        bin_center_1 = cluster_bin_map[i]
                    else:
                        if ndim == 1:
                            bin_center_1 = bin_centers[bin_map[i]]
                        if ndim == 2:
                            bin_center_1 = bin_centers[bin_map[i]][0]
                            bin_center_2 = bin_centers[bin_map[i]][1]
                    if cluster_binning:
                        fep_pb_output_file.write("  %12d %12d " % (bin_map[i], bin_center_1))
                    else:
                        if ndim == 1:
                            fep_pb_output_file.write("  %12d %12.3f " % (bin_map[i], bin_center_1))
                        if ndim == 2:
                            fep_pb_output_file.write("  %12d %12.3f %12.3f " % (bin_map[i], bin_center_1, bin_center_2))
                    # if bin clustering, include pmf information in feppb file
                    if cluster_binning:
                        fep_pb_output_file.write("%12.3f %12.3f %12.3f %12.3f " % (f_i[i], df_i[i], E_i[i], S_i[i]))

                    for j in range(len(reduced_fep_data)):
                        fep_pb_output_file.write("%12.3f " % fep_values[i][j])
                    fep_pb_output_file.write("\n")
                fep_pb_output_file.close()

        # Write out PMF for reduced bins.
        pmf_output_file.write("# PMF, energy and entropy for reduced (non-empty) bins (in units of kT)\n")
        if ndim == 1:
          pmf_output_file.write("# %12s %12s %12s %12s %12s %12s\n" % ('bin', 'bin_center_1', 'f', 'df', 'e', 's'))
        elif ndim == 2:
          pmf_output_file.write("# %12s %12s %12s %12s %12s %12s %12s\n" % ('bin', 'bin_center_1', 'bin_center_2', 'f', 'df', 'e', 's'))
        for i in range(nreduced_bins):
            if cluster_binning:
                bin_center_1 = cluster_bin_map[i]
            else:
                if ndim == 1:
                    bin_center_1 = bin_centers[bin_map[i]]
                if ndim == 2:
                    bin_center_1 = bin_centers[bin_map[i]][0]
                    bin_center_2 = bin_centers[bin_map[i]][1]
            if cluster_binning:
                pmf_output_file.write("  %12d %12d %12.3f %12.3f %12.3f %12.3f\n" % (bin_map[i], bin_center_1, f_i[i], df_i[i], E_i[i], S_i[i]))
            else:
                if ndim == 1:
                    pmf_output_file.write("  %12d %12.3f %12.3f %12.3f %12.3f %12.3f\n" % (bin_map[i], bin_center_1, f_i[i], df_i[i], E_i[i], S_i[i]))
                if ndim == 2:
                    pmf_output_file.write("  %12d %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f\n" % (bin_map[i], bin_center_1, bin_center_2, f_i[i], df_i[i], E_i[i], S_i[i]))

        pmf_output_file.close()

    ev_output_file.close()
    cv_output_file.close()
