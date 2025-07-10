'''
Author: <Logan Vossler>

==========================================================
Make Site Frequency Spectra
==========================================================

Description:
This script uses dadi to make and export 1D or 2D Allele Frequency Spectra files/plots.
It will also save bootstrapped SFS data for later use in dadi uncertainty analysis.

File Requirements:
    - Base Parameter File from Genomics Main Lab Repository
    - Dadi-Specific Parameter File from Genomics Main Lab Repository
    - A previously generated data dictionary .pkl file present in result directory
    - POP ID File specifying the population id for each individual in the VCF (See dadi-documentation for info on formatting)
'''


# Required Modules
#==========================================================
import dadi, random, os
import matplotlib.pyplot as plt
import dill as pkl
import argparse


# CMD-LINE Arguments
#==========================================================
parser=argparse.ArgumentParser()
parser.add_argument("-f", "--folder_name", type=str, help='The name of the folder (under your general OUTDIR) you want to generate containing your results')
parser.add_argument("-p","--pop_ids", type=str, help='List of Pop_IDs Params', nargs='+')
parser.add_argument("-n", "--num_chroms", type=str, help='List of ints representing number of chromosomes per population', nargs='+')
parser.add_argument("-l", "--lowpass", type=bool, help='Do Lowpass Pipeline')
parser.add_argument("-o","--outdir", type=str, help='Out Directory Path')
parser.add_argument("-b","--bootparams", type=str, help='List of bootstrapping params', nargs='+', default="100 1e-7")
parser.add_argument("-t","--polarize", type=bool, help='State this parameter to UNFOLD the SFS')
args = parser.parse_args()


# Function Definitions
#==========================================================
def bootstrap(dd, pop_ids, num_chrom, result_dir, Nboot=100, chunk_size=1e7):
    '''
    Create bootstrapped datasets from our SNP data dictionary.
    Parameters:
        dd: A Data Dictionary
        pop_ids: List of strings of spp/pop names for which we are generating bootstrapped datasets
        num_chrom: List of ints representing the number of chromosomes for each species we are bootstrapping
        result_dir: String representing dadi results directory path to where the bootstraps will be saved
        Nboot: Int representing the number of bootstrapped SFS to create (defaults to 100)
        chunk_size: Int (can be scientific notation) representing the chunk size of each bootstrapped SFS (defaults to 1e7)
    Returns:
        None
    '''
    # Break data dictionary into chunks (list of dictionary genome chunks)
    chunks = dadi.Misc.fragment_data_dict(dd, chunk_size)
    # Set any random seed so that the same random genome chunks are selected for non/synonymous mutations
    random.seed(1762)
    # Get a list of sfs from bootstrapped genomes for each pop combo
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
    Takes an sfs object and constructs a plot (either 1D or 2D).
    If the plot is 2D, Pop0 will plot on the yaxis and Pop1 will plot on the xaxis. 
    Parameters:
        sfs: dadi SFS object
        result_dir: String representing the dadi results directory path to where plot will be saved
        pop_ids: pop_ids: List of strings of pop names
        dim: Int representing dimensionality of the sfs plot being created (defaults to 2)
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
    Saves a .pkl file of a depth-of-coverage distribution for a given population.
    This coverage distribution is only needed if using LowPass for Low Coverage data in future model-making.
    We generate the file here (instead of in dadi_2_demo_model.py) to avoid having to load in the large data dictionary multiple times (thus greatly reducing computational resources)
    Parameters:
        dd: A Data Dictionary
        sfs: dadi SFS object
        pop_ids: List of strings of pop names
        dim: Int representing the dimensionality of the sfs
    '''
    cov_dist = LowPass.compute_cov_dist(dd, pop_ids)
    with open(lowpass_dir + '_'.join(pop_ids) + '_cov_dist.pkl', 'wb') as file:
        pkl.dump(cov_dist, file)


# Main
#==========================================================
'''
1) Clean some CMD-Line Parameters
2) Verify/Make Result Directory
3) Load Data Dictionary
4) Generate, Plot, and Save SFS
5) Generate Cov-Dist if needed
6) Generate Bootstrapped SFS
'''
#========================================
# Clean up some arguments
pop_ids = args.pop_ids[0].split()
num_boots = int(args.bootparams[0].split()[0])
chunk_size = float(args.bootparams[0].split()[1])
num_chroms = [int(num) for num in args.num_chroms[0].split()]

#========================================
# Check if dadi-specific results directory exists in specified outdir. If not, create it.
print('Verifying Directories...')
# We will store our SFS data in a general results folder for each population combo...
result_dir = args.outdir + args.folder_name + '/'
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

#...but if using lowpass in models later, also make a lowpass directory inside specified results folder
if args.lowpass:
    from dadi.LowPass import LowPass
    print('Generating Lowpass directory...')
    lowpass_dir = result_dir + '/lowpass/'
    if not os.path.exists(lowpass_dir):
        os.makedirs(lowpass_dir)

#========================================
# Make Data dictionary and save to file.
if os.path.exists(result_dir + 'dd.pkl'):
    print('Data Dictionary .pkl file found in result directory. Loading this file into job...')
    with open(result_dir + 'dd.pkl', 'rb') as file:
        dd = pkl.load(file)
else:
    print('Error: No data dictionary found in directory')

#========================================
# Make Spectrum objects
print('Generating SFS for ' + ' '.join(pop_ids) + '...')
data_fs = dadi.Spectrum.from_data_dict(dd, pop_ids, polarized=args.polarize, projections=num_chroms)

# A Hack to determine fs dimensionality (If first element in fs array is NOT an array/spectrum object, then it is 1D)
# Thought there was a native dadi method to get dimensionality, but I guess not...
dim = 2 if type(data_fs[0]) == dadi.Spectrum_mod.Spectrum else 1

# Plot the SFS and save plot to file
print('Plotting SFS for ' + ' '.join(pop_ids) + '...')
plot_sfs(data_fs, result_dir, pop_ids, dim)

# Save SFS to files
print('Saving SFS for ' + ' '.join(pop_ids) + '...')
data_fs.to_file(result_dir + '_'.join(pop_ids) + '_fs')

# If doing lowpass workflow, save coverage distribution for future modeling
if args.lowpass:
    print('Saving lowpass coverage distribution to temp file in lowpass dir...')
    save_cov_dist(dd, lowpass_dir, pop_ids, dim)

# Make Bootstrapped SFS files
print('Generating Bootstrapped SFS for ' + ' '.join(pop_ids) + '...')
bootstrap(dd, pop_ids, num_chroms, result_dir, num_boots, chunk_size)

print('\n**SFS Creation Complete**')

