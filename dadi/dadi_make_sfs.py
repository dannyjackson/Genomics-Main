'''
Author: <Logan Vossler>

==========================================================
Make 2-Dimensional Site Frequency Spectra
==========================================================

Description:
This script uses dadi to make and export 2-Dimensional Allele Frequency Spectrum files and plots.
It will also save bootstrapped SFS data for later use in dadi uncertainty analysis.

File Requirements:
    - Base Parameter File from Genomics Main Lab Repository
    - Dadi-Specific Parameter File from Genomics Main Lab Repository
    - VCF File containing SNP data for desired populations for analysis
    - POP ID File specifying the population id for each individual in the VCF (See dadi-documentation for info on formatting)
'''


# Required Modules
#==========================================================
import dadi, random, os, sys, json5
import matplotlib.pyplot as plt
import pickle as pkl
from pathlib import Path

# Function Definitions
#==========================================================
def bootstrap(dd, pop_ids, num_chrom, result_dir):
    '''
    This function creates bootstrapped datasets from our SNP data dictionary.
    Parameters:
        dd: A Data Dictionary
        pop_ids: A 2 element list containing strings of species names for which we are generating bootstrapped datasets
        num_chrom: A 2 element list containing integers representing the number of chromosomes for each species we are bootstrapping
        result_dir: A string representing the dadi results directory path to where the bootstraps will be saved
    Returns:
        None
    '''
    # State number of bootstrapped datasets desired and genome chunks
    Nboot, chunk_size = 100, 1e7
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

def plot_sfs(sfs, pop_ids, result_dir):
    '''
    This function takes an 2D sfs object and constructs a plot to be saved to the dadi_results directory. 
    Pop0 will plot on the yaxis and Pop1 will plot on the xaxis. 
    Parameters:
        sfs: A 2D sfs object
	    pop_ids: A 2 element list containing strings for species being plotted. ###yaxis is the first list element. xaxis is the second.####
        result_dir: A string representing the dadi results directory path to where the bootstraps will be saved
    Returns:
        None
    '''
    plot_spectrum = dadi.Plotting.plot_single_2d_sfs(sfs, vmin=1, pop_ids=pop_ids)
    plt.savefig(result_dir + '_'.join(pop_ids) + '_2d_spectrum.png')
    plt.clf()


def main():
    '''
    1) Import Base Parameters
    2) Import Dadi-SFS Specific Parameters
    3) Check for required directories
    4) Make Data Dictionary
    '''
    #========================================
    # Import base parameters from user-inputted params_base.sh file
    base_params = Path(sys.argv[1]).read_text().strip().split('\n')
    for line in base_params:
        if 'OUTDIR=' in line: outdir = line.split('=')[1].split('#')[0].strip()
        if 'PROGDIR=' in line: progdir = line.split('=')[1].split('#')[0].strip()
        if 'PROJHUB=' in line: github = line.split('=')[1].split('#')[0].strip()
        if 'SCRIPTDIR=' in line: scriptdir = line.split('=')[1].split('#')[0].strip()

    # Import dadi-specific parameters for sfs creation from user-inputted dadi_params.json file
    with open(sys.argv[2], 'r') as file:
        dadi_params = json5.load(file)
    vcffile = dadi_params['VCF PATH']
    popfile = dadi_params['POP PATH']
    sfsparams = dadi_params['SFS PARAMS']
    #========================================
    # Check if dadi-specific results directory exists in specified outdir. If not, create it.
    result_dir = outdir + 'dadi_results/'
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    # Make Data dictionary and save to file.
    dd = dadi.Misc.make_data_dict_vcf(vcffile, popfile, calc_coverage=True)
    with open('dadi_results/dd.bpkl', 'wb') as file:
        pkl.dump(dd, file, 2)

    # Enter a loop to perform tasks below for each species combo
    for dct in sfsparams:
        pop_ids = sfsparams[dct][0]
        num_chrom = sfsparams[dct][1]
        polarize = sfsparams[dct][2]

        # Make Spectrum objects
        data_fs = dadi.Spectrum.from_data_dict(dd, pop_ids, polarized=polarize, projections=num_chrom)
        
        # Plot the SFS and save plot to file
        plot_sfs(data_fs, pop_ids, result_dir)
        
        # Save SFS to files
        data_fs.to_file(result_dir + '_'.join(pop_ids) + '_fs')
        
        # Make Bootstrapped SFS files
        bootstrap(dd, pop_ids, num_chrom, result_dir)


if __name__ == '__main__':
    main()
