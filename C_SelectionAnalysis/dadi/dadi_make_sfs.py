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
    #boots = dadi.Misc.bootstraps_from_dd_chunks(chunks, Nboot, pop_ids=pop_ids, polarized=False, projections=num_chrom)
    boots = []
    for i in range(Nboot):
        chosen = random.choice(chunks, k=len(chunks))
        temp = {}
        for j in chosen:
            temp.update(j)
        boots.append(dadi.Spectrum.from_data_dict(temp, pop_ids, num_chrom))

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
    print('\nMaking Data Dictionary...')
    dd = dadi.Misc.make_data_dict_vcf(vcffile, popfile, calc_coverage=True)
    print('Saving Data Dictionary to results directory...')
    with open(result_dir + 'dd.pkl', 'wb') as file:
        pkl.dump(dd, file)
    
    #========================================
    # Enter a loop to make SFS for each species combo
    #for dct in sfsparams:
        #print('\nGetting ' + dct + ' Parameters...')
        #pop_ids, num_chrom, polarize = sfsparams[dct]

        # Make Spectrum objects
        #print('Generating SFS for ' + dct + '...')
        #data_fs = dadi.Spectrum.from_data_dict(dd, pop_ids, polarized=polarize, projections=num_chrom)
        
        # Plot the SFS and save plot to file
        #print('Plotting SFS for ' + dct + '...')
        #plot_sfs(data_fs, pop_ids, result_dir)
        
        # Save SFS to files
        #print('Saving SFS for ' + dct + '...')
        #data_fs.to_file(result_dir + '_'.join(pop_ids) + '_fs')
        
        # Make Bootstrapped SFS files
        #print('Generating Bootstrapped SFS for ' + dct + '...')
        #bootstrap(dd, pop_ids, num_chrom, result_dir, num_boots, chunk_size)

    print('\n**SFS Creation Complete**')


if __name__ == '__main__':
    main()
