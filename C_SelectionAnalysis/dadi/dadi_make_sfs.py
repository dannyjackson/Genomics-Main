'''
Author: <Logan Vossler>

==========================================================
Make Site Frequency Spectra
==========================================================

Description:
This script uses dadi to make and export 1 or 2 Dimensional Allele Frequency Spectrum files and plots.
It will also save bootstrapped SFS data for later use in dadi uncertainty analysis.

File Requirements:
    - Base Parameter File from Genomics Main Lab Repository
    - Dadi-Specific Parameter File from Genomics Main Lab Repository
    - VCF File containing SNP data for desired populations for analysis
    - POP ID File specifying the population id for each individual in the VCF (See dadi-documentation for info on formatting)
'''


# Required Modules
#==========================================================
import dadi, random, os, sys
import matplotlib.pyplot as plt
import dill as pkl
from dadi.LowPass import LowPass # Can comment out if not using LowPass workflow
from pathlib import Path
import json5 # Can switch to normal json module if this one causes issues

# Function Definitions
#==========================================================
def bootstrap(dd, pop_ids, num_chrom, result_dir, Nboot=100, chunk_size=1e7):
    '''
    This function creates bootstrapped datasets from our SNP data dictionary.
    Parameters:
        dd: A Data Dictionary
        pop_ids: A 2 element list containing strings of species names for which we are generating bootstrapped datasets
        num_chrom: A 2 element list containing integers representing the number of chromosomes for each species we are bootstrapping
        result_dir: A string representing the dadi results directory path to where the bootstraps will be saved
        Nboot: An integer representing the number of bootstrapped SFS to create (defaults to 100)
        chunk_size: An integer (can be scientific notation) representing the chunk size of each bootstrapped SFS (defaults to 1e7)
    Returns:
        None
    '''
    # Break data dictionary into chunks (list of dictionary genome chunks)
    chunks = dadi.Misc.fragment_data_dict(dd, chunk_size)
    # Set any random seed so that the same random genome chunks are selected for non/synonymous mutations
    random.seed(1762)
    # Get a list containing sfs from bootstrapped genomes for each pop combo
    boots = dadi.Misc.bootstraps_from_dd_chunks(chunks, Nboot, pop_ids=pop_ids, polarized=False, projections=num_chrom)

    # Check for Directories for dadi sfs bootstraps
    if not os.path.exists(result_dir + 'bootstraps/'):
        os.makedirs(result_dir + 'bootstraps/')
    if not os.path.exists(result_dir + 'bootstraps/' + '_'.join(pop_ids) + '/'):
        os.makedirs(result_dir + 'bootstraps/' + '_'.join(pop_ids) + '/')

    # Save bootstrapped SFS to files
    boot_dir = result_dir + 'bootstraps/' + '_'.join(pop_ids) + '/'
    for i in range(len(boots)):
        boots[i].to_file(boot_dir + '_'.join(pop_ids) + 'boots{0}.fs'.format(str(i)))

def plot_sfs(sfs, result_dir, pop_ids, dim=2):
    '''
    This function takes an sfs object and constructs a plot to be saved to the dadi_results directory. 
    It can create 1 and 2 dimensional sfs
    If the plot is 2D, Pop0 will plot on the yaxis and Pop1 will plot on the xaxis. 
    Parameters:
        sfs: A 2D sfs object
        result_dir: A string representing the dadi results directory path to where the bootstraps will be saved
        pop_ids: pop_ids: A list containing strings of pop names
        dim: An integer representing the dimensionality of the sfs plot being created
    Returns:
        None
    '''
    if dim == 1:
        # Plot 1d SFS
        spectrum_plot = dadi.Plotting.plot_1d_fs(sfs)
        plt.savefig(result_dir + pop_ids[0] + '_1d_spectrum.png')
    elif dim == 2:
        # Plot 2d SFS
        spectrum_plot = dadi.Plotting.plot_single_2d_sfs(sfs, vmin=1, pop_ids=pop_ids)
        plt.savefig(result_dir + '_'.join(pop_ids) + '_2d_spectrum.png')
    plt.clf()

def save_cov_dist(dd, result_dir, pop_ids, dim=2):
    '''
    This function saves a .pkl file of a depth-of-coverage distribution for a given population.
    This coverage distribution is only needed if using LowPass for Low Coverage data in future model-making.
    We generate the file here (instead of in dadi_make_2dmodel.py) to avoid having to load in the large data dictionary
    multiple times (thus greatly reducing computational resources)
    Parameters:
        dd: A Data Dictionary
        sfs: A 2D SFS object
        pop_ids: pop_ids: pop_ids: A list containing strings of pop names
        dim: An integer representing the dimensionality of the sfs plot being created
    '''
    fname = pop_ids[0] if dim == 1 else '_'.join(pop_ids)
    cov_dist = LowPass.compute_cov_dist(dd, pop_ids)
    with open(result_dir + fname + '_cov_dist.pkl', 'wb') as file:
        pkl.dump(cov_dist, file)


# Main
#==========================================================
def main():
    '''
    1) Import Base Parameters
    2) Import Dadi-SFS Specific Parameters
    3) Check for required directories
    4) Make Data Dictionary
    5) Generate and Plot SFS
    6) Generate Bootstrapped SFS
    '''
    #========================================
    # Basic check for enough parameter files inputted
    print('Checking File Arguments...')
    if len(sys.argv) != 3:
        raise IndexError('Not enough arguments. Did you include both the Base and dadi-specific parameter files?')
    
    # Import base parameters from user-inputted params_base.sh file
    print('Storing Needed Base Parameters...')
    base_params = Path(sys.argv[1]).read_text().strip().split('\n')
    outdir = None
    for line in base_params:
        if 'OUTDIR=' in line: outdir = line.split('=')[1].split('#')[0].strip() 
    if not outdir:
        raise NameError("Couldn't find Out-Directory. Check that OUTDIR is specified as in example Genomics-Main base_params.sh.")
    
    # Import dadi-specific parameters for sfs creation from user-inputted dadi_params.json file
    print('Storing Needed dadi Parameters...')
    with open(sys.argv[2], 'r') as file:
        dadi_params = json5.load(file)
    job_name = dadi_params['JOB NAME']
    vcffile = dadi_params['VCF PATH']
    popfile = dadi_params['POP PATH']
    sfsparams = dadi_params['SFS PARAMS']
    num_boots, chunk_size = dadi_params['BOOTSTRAP PARAMS']
    lowpass = dadi_params['LOWPASS']
    
    #========================================
    # Check if dadi-specific results directory exists in specified outdir. If not, create it.
    print('Verifying Directories...')
    # If using lowpass, make a lowpass directory inside specified results folder
    result_dir = outdir + job_name + '/lowpass/' if lowpass else outdir + job_name + '/'
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    #========================================
    # Make Data dictionary and save to file.
    if os.path.exists(result_dir + 'dd.pkl'):
        print('Data Dictionary .pkl file found in result directory. Loading this file into job...')
        with open(result_dir + 'dd.pkl', 'rb') as file:
            dd = pkl.load(file)
    else:
        print('\nMaking Data Dictionary...')
        dd = dadi.Misc.make_data_dict_vcf(vcffile, popfile, calc_coverage=True)
        print('Saving Data Dictionary to results directory...')
        with open(result_dir + 'dd.pkl', 'wb') as file:
            pkl.dump(dd, file)
    #========================================
    # Enter a loop to make SFS for each species combo
    for dct in sfsparams:
        print('\nGetting ' + dct + ' Parameters...')
        pop_ids, num_chrom, polarize, dim = sfsparams[dct]

        # Make Spectrum objects
        print('Generating SFS for ' + dct + '...')
        data_fs = dadi.Spectrum.from_data_dict(dd, pop_ids, polarized=polarize, projections=num_chrom)
        
        # Plot the SFS and save plot to file
        print('Plotting SFS for ' + dct + '...')
        plot_sfs(data_fs, result_dir, pop_ids, dim)
        
        # If doing lowpass workflow, save coverage distribution for future modeling
        if lowpass:
            print('Saving lowpass coverage distribution to temp file...')
            save_cov_dist(dd, data_fs, result_dir, pop_ids, dim)
        
        # Save SFS to files
        print('Saving SFS for ' + dct + '...')
        data_fs.to_file(result_dir + '_'.join(pop_ids) + '_fs')
        
        # Make Bootstrapped SFS files
        print('Generating Bootstrapped SFS for ' + dct + '...')
        bootstrap(dd, pop_ids, num_chrom, result_dir, num_boots, chunk_size)

    print('\n**SFS Creation Complete**')


if __name__ == '__main__':
    main()
