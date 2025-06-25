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
import argparse

print('define arguments')
parser=argparse.ArgumentParser()
parser.add_argument("-j", "--job_name", type=str, help='Name of Job')
parser.add_argument("-f", "--folder_name", type=str, help='The name of the folder (under your general OUTDIR) you want to generate containing your results')
parser.add_argument("-p","--pop_ids", type=str, help='List of Pop_IDs Params', nargs='+')
parser.add_argument("-n", "--num_chroms", type=str, help='List of ints representing number of chromosomes per population', nargs'+')
parser.add_argument("-l", "--lowpass", type=bool, help='Do Lowpass Pipeline')
parser.add_argument("-o","--outdir", type=str, help='Out Directory Path')
parser.add_argument("-v","--vcfpath", type=str, help='Path to vcf file')
parser.add_argument("-i","--poppath", type=str, help='Path to pop assignment file')
parser.add_argument("-b","--bootparams", type=str, help='List of bootstrapping params', nargs='+', default="100 1e-7")
parser.add_argument("-t","--polarize", type=bool, help='State this parameter to UNFOLD the SFS')
args = parser.parse_args()

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

def save_cov_dist(dd, lowpass_dir, pop_ids, dim=2):
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
    with open(lowpass_dir + fname + '_cov_dist.pkl', 'wb') as file:
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
    # Clean up some arguments
    pop_ids = args.pop_ids[0].split()
    num_boots, chunk_size = args.bootparams[0].split()
    num_chroms = [int(num) for num in args.num_chroms[0].split()]
    print(pop_ids)
    print(num_boots, chunk_size)
    print(num_chroms)
    
    #========================================
    # Check if dadi-specific results directory exists in specified outdir. If not, create it.
    print('Verifying Directories...')
    # We will store our SFS data in a general results folder for each population combo...
    result_dir = args.outdir + args.folder_name + '/'
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    #...but if using lowpass in models later, also make a lowpass directory inside specified results folder
    if args.lowpass:
        print('Generating Lowpass directory...')
        lowpass_dir = result_dir + '/lowpass/'
        if not os.path.exists(lowpass_dir):
            os.makedirs(lowpass)

    #========================================
    # Make Data dictionary and save to file.
    if os.path.exists(result_dir + 'dd.pkl'):
        print('Data Dictionary .pkl file found in result directory. Loading this file into job...')
        with open(result_dir + 'dd.pkl', 'rb') as file:
            dd = pkl.load(file)
    else:
        print('\nMaking Data Dictionary...')
        # First, if not doing lowpass, then don't bother including coverage info in dd
        cov = True if args.lowpass else False
        dd = dadi.Misc.make_data_dict_vcf(args.vcfpath, args.poppath, calc_coverage=cov)
        print('Saving Data Dictionary to results directory...')
        with open(result_dir + 'dd.pkl', 'wb') as file:
            pkl.dump(dd, file)
    #========================================
    # Make Spectrum objects
    print('Generating SFS for ' + ' '.join(pop_ids) + '...')
    data_fs = dadi.Spectrum.from_data_dict(dd, pop_ids, polarized=args.fold, projections=num_chroms)
    
    # Plot the SFS and save plot to file
    print('Plotting SFS for ' + ' '.join(pop_ids) + '...')
    plot_sfs(data_fs, result_dir, pop_ids)
    
    # If doing lowpass workflow, save coverage distribution for future modeling
    if args.lowpass:
        print('Saving lowpass coverage distribution to temp file in lowpass dir...')
        save_cov_dist(dd, lowpass_dir, pop_ids)
    
    # Save SFS to files
    print('Saving SFS for ' + dct + '...')
    data_fs.to_file(result_dir + '_'.join(pop_ids) + '_fs')
    
    # Make Bootstrapped SFS files
    print('Generating Bootstrapped SFS for ' + dct + '...')
    bootstrap(dd, pop_ids, num_chroms, result_dir, num_boots, chunk_size)

    print('\n**SFS Creation Complete**')


if __name__ == '__main__':
    main()
