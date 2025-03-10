In terms of installing Timesweeper, I was able to do this using a python virtual environment (information on these can be found here: https://hpcdocs.hpc.arizona.edu/software/popular_software/python/#installing-python-packages-using-a-virtual-environment). For example:

module load python/3.8/3.8.2
python3 -m venv --system-site-packages timesweeper_env

source timesweeper_env/bin/activate
module load slim samtools bcftools 
git clone https://github.com/SchriderLab/Timesweeper.git
cd Timesweeper/
# edit setup.py to be:
from setuptools import setup, find_packages

if __name__ == "__main__":
    with open("README.md", "r") as readme:
        readme_txt = readme.read()
    
    setup(name="timesweeper",packages=find_packages(),long_description=readme_txt, long_description_content_type='text/markdown')

pip install --no-cache-dir .
pip install -r requirements.txt --no-cache-dir 
pip list --format=freeze > /xdisk/mcnew/finches/dannyjackson/timesweeper_test/updated_requirements.txt

pip install scikit-allel==1.3.5 --no-cache-dir 
# Timesweeper test run

echo 'export PATH="/home/u15/dannyjackson/.local/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

cd /xdisk/mcnew/finches/dannyjackson/timesweeper_test/examplerun

interactive -a mcnew -t 6:00:00 -g

# module load slim samtools bcftools cuda tensorflow/nvidia/2.9.1


# custom simulation script
timesweeper sim_custom -y example_config.yaml 

# compute average missingness of samples
awk '{ total += $6 } END { print total/NR }' /xdisk/mcnew/finches/dannyjackson/finches/summarystats/plink.imiss
0.0218649

# Make Training Data (condense)
timesweeper condense -o test.pkl -m 0.02 -y  example_config.yaml 

# Neural networks
cd /xdisk/mcnew/finches/dannyjackson/timesweeper_test/programs/Timesweeper/build/lib/timesweeper

awk '800<NR && NR<825' /home/u15/dannyjackson/.local/lib/python3.8/site-packages/tensorflow/python/trackable/data_structures.py

timesweeper train -i test.pkl -y example_config.yaml 
# prepare vcf
# detect sweeps

sed '/##INFO=<ID=VDB/s/,Version="3"//' /xdisk/mcnew/finches/ljvossler/finches/test_genom_main_dadi_scripts/one_chrom_qualitysort.vcf > test.vcf

timesweeper detect -i test.vcf -y example_config.yaml -o /xdisk/mcnew/finches/dannyjackson/timesweeper_test/examplerun/output
