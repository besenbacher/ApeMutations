The workflow is created using gwf:
http://gwf.readthedocs.io/en/latest/

gwf can be installed using conda (see link above) or alternativly from source:
https://github.com/gwforg/gwf

Besides gwf samtools and bwa also needs to be installed to run the workflow. And
they should be in the PATH of the user running the jobs (otherwise the workflow
should be updated with the full path to the programs).

The path to pairs of fastq files should be specified in the file called
fastq_files.txt. There can be multiple fastq pairs for each individual.

The reference genome or (a dynamic link to it) should be placed in the subfolder
called ref-genomes. It can be downloaded here:
http://hgdownload.soe.ucsc.edu/goldenPath/panTro5/bigZips/panTro5.fa.gz

Once gwf is installed you should run the following command in this directory to
tell it that jobs should be submitted using slurm:

gwf config set backend slurm

To submit all jobs just type:

gwf run

To check status of jobs type:

gwf status

To check error message of finished job type:

gwf -e logs {job_name}
