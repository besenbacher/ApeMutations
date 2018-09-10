######################################################################################
#A.SET ENVIRONMENT AND IMPORT SYS,OS
###################################################################################### 
from gwf import Workflow
import sys, os
import itertools
from itertools import *
gwf = Workflow()

######################################################################################
#B.FUNCTIONS AND TEMPLATES
######################################################################################

# Upstream analysis: map reads and merge bams for individual

def get_index_bwa(ref, ref_genome_fasta, ref_genome_am, ref_genome_ann, ref_genome_pac, ref_genome_sa, completed):
  inputs = [ref_genome_fasta]
  outputs = [ref_genome_am, ref_genome_ann, ref_genome_pac, ref_genome_sa, completed]
  options = {"memory": "8g","walltime":"12:00:00"}
  spec='''
  source activate prj1
  source /com/extra/bwa/LATEST/load.sh

  bwa index -p {ref} -a bwtsw {ref_genome_fasta} 

  echo "Completed at "$(date) > {completed}'''.format(ref = ref, ref_genome_fasta = ref_genome_fasta, completed = completed)

  return inputs, outputs, options, spec

def map_bwa(ref, read_1, read_2, ref_genome_am, ref_genome_ann, ref_genome_pac, ref_genome_sa, bam, completed):
  inputs = [read_1, read_2, ref_genome_am, ref_genome_ann, ref_genome_pac, ref_genome_sa]
  outputs = [bam, completed]
  options = {"memory": "64g","walltime":"100:00:00","cores":16}
  spec='''
  source activate prj1
  source /com/extra/bwa/LATEST/load.sh

  bwa mem -t 16 {ref} {read_1} {read_2} | samtools sort | samtools rmdup -s - {bam} 

  echo "Completed at "$(date) > {completed}'''.format(ref = ref, read_1 = read_1, read_2 = read_2, bam=bam, completed = completed)

  return inputs, outputs, options, spec

def merge_bams(inputbams, merged_bam, merged_bam_rg, name, prev_completed, completed):
  inputs = [prev_completed]
  outputs = [merged_bam, merged_bam_rg, completed]
  options = {"memory": "8g","walltime":"100:00:00"}
  spec='''
  source activate prj1
  source /com/extra/bwa/LATEST/load.sh

  rg_file=merged_bam_files/{name}.rg.txt
  rm -f $rg_file
  for input in {inputbams}; do
    bn=`basename $input .bam`
    echo -e "@RG\tID:$bn\tSM:{name}\tLB:{name}\tPL:Illumina" >> $rg_file
  done
  samtools merge -rh $rg_file - {inputbams} | samtools rmdup -s - {merged_bam}

  echo "Completed at "$(date) > {completed}'''.format(inputbams = inputbams, name=name, merged_bam = merged_bam, completed=completed)
  
  return inputs, outputs, options, spec

def read_family_description(fname):
    f = open(fname)
    res = []
    for line in f:
        L = line.split()
        res.append((L[0],L[1],L[2]))
    f.close()
    return res

### Helper function to add extra 'input' and 'output' files to a options dict
def add_options(options, extra):
    for inout in ['input', 'output']:
        if inout in extra:
            if inout in options:
                if type(options[inout]) != list:
                    options[inout] = [options[inout]]
                if type(extra[inout]) != list:
                    extra[inout] = [extra[inout]]
                options[inout] = list(set(options[inout] + extra[inout]))
            else:
                options[inout] = extra[inout]
    return options

def picard_index_reference(ref_genome_fasta,ref_genome_dict,ref_genome_fai,completed):
  inputs = [ref_genome_fasta]
  outputs = [ref_genome_dict,ref_genome_fai,completed]
  options = {"memory": "1g","walltime":"01:00:00"}
  spec='''
  source activate prj1
  source /com/extra/samtools/LATEST/load.sh
  
  samtools faidx ref-genomes/{ref_genome_fasta}

  source /com/extra/java/8/load.sh
  source /com/extra/picard/LATEST/load.sh
  
  picard CreateSequenceDictionary R=ref-genomes/{ref_genome_fasta} O=ref-genomes/{ref_genome_dict}

  echo "Completed at "$(date) > {completed}'''.format(ref_genome_fasta=ref_genome_fasta,ref_genome_dict=ref_genome_dict, completed=completed)
  
  return inputs, outputs, options, spec

def samtools_index_bam(bam,bai,completed):
  inputs = [bam]
  outputs = [bai,completed]
  options = {"memory": "5g","walltime":"10:00:00"}
  spec='''
  source activate prj1
  source /com/extra/samtools/LATEST/load.sh
  
  samtools index {bam}

  echo "Completed at "$(date) > {completed}'''.format(bam=bam, completed=completed)
  
  return inputs, outputs, options, spec

# Mark duoplicates ? 
      #java -jar picard.jar MarkDuplicates \
      #I=input.bam \
      #O=marked_duplicates.bam \
      #M=marked_dup_metrics.txt

# before/after plots to visualize the effects of the recalibration process
# if we don't have high quality vcf files for Great Apes:
#bootstrap a set of known variants: 
  #First do an initial round of variant calling on your original, unrecalibrated data.
  #Then take the variants that you have the highest confidence in and use that set as the database of known variants by feeding it as a VCF file to the BaseRecalibrator.
  #Finally, do a real round of variant calling with the recalibrated data. These steps could be repeated several times until convergence.

# BaseRecalibrator: more than one vcf can be used

def calc_recalibrate_info(ref_genome_fasta, ref_genome_dict, merged_bam, merged_bai, knownlist, recalibration_report, completed):
  inputs = [ref_genome_fasta, ref_genome_dict, merged_bam, merged_bai, knownlist]
  outputs = [recalibration_report, completed]
  options = {"cores": 16}
  spec='''
  source activate prj1
  source /com/extra/java/8/load.sh
  source /com/extra/GATK/4/load.sh

 gatk BaseRecalibrator \
  -nct 16 \
  -I {merged_bam} \
  -R {ref_genome_fasta} \
  --known-sites {knownlist} \
  -O {recalibration_report}
  
  echo "Completed at "$(date) > {completed}'''.format(ref_genome_fasta=ref_genome_fasta, merged_bam = merged_bam, knownlist=knownlist, recalibration_report=recalibration_report, completed=completed)
  
  return inputs, outputs, options, spec

# change temporal dir path

def recalibrate(ref_genome_fasta, ref_genome_dict, merged_bam, merged_bai, recalibration_report, recalibrated_bam, recalibrated_bai, completed):
  inputs = [ref_genome_fasta, ref_genome_dict, merged_bam, merged_bai, recalibration_report]
  outputs = [recalibrated_bam, recalibrated_bai, completed]
  options = {"cores": 16, "memory":"50g"}
  spec='''
  source activate prj1
  source /com/extra/java/8/load.sh
  source /com/extra/GATK/4/load.sh
  
  java -Djava.io.tmpdir=/scratch/$PBS_JOBID/tmp \
    -Djava.awt.headless=true \
    -Xmx50g \
    -jar GenomeAnalysisTK.jar \
    -T PrintReads \
    -R {ref_genome_fasta} \
    -I {merged_bam} \
    -BQSR {recalibration_report}\
    -o {recalibrated_bam}

  echo "Completed at "$(date) > {completed}'''.format(ref_genome_fasta=ref_genome_fasta, ref_genome_dict=ref_genome_dict, merged_bam=merged_bam, merged_bai=merged_bai, recalibration_report=recalibration_report, recalibrated_bam=recalibrated_bam, recalibrated_bai=recalibrated_bai, completed=completed)

  return inputs, outputs, options, spec

# Same read filters?

def filter_bam_files(recalibrated_bam, recalibrated_bai, ref_genome_fasta, filtered_bam, completed):
  inputs = [recalibrated_bam, recalibrated_bai, ref_genome_fasta]
  outputs = [filtered_bam, completed]
  options = {"cores": 16, "walltime":"12:00:00"}
  spec='''
  source activate prj1
  source /com/extra/java/8/load.sh
  source /com/extra/GATK/4/load.sh

  gatk \
     -R  {ref_genome_fasta} \
     -T PrintReads \
     -I {recalibrated_bam} \
     -o {filtered_bam} \
     -nct 16 \
     --read_filter BadCigar \
     --read_filter DuplicateRead \
     --read_filter FailsVendorQualityCheck \
     --read_filter HCMappingQuality \
     --read_filter MappingQualityUnavailable \
     --read_filter NotPrimaryAlignment \
     --read_filter UnmappedRead \
     --filter_bases_not_stored \

  echo "Completed at "$(date) > {completed}'''.format(recalibrated_bam=recalibrated_bam, recalibrated_bai=recalibrated_bai, ref_genome_fasta=ref_genome_fasta, filtered_bam=filtered_bam, completed=completed)

  return inputs, outputs, options, spec

# Maybe have to change to:
   #gatk --java-options "-Xmx4g" HaplotypeCaller  \
   #-R {ref_genome_fasta}  \
   #-I input.bam \
   #-O output.g.vcf.gz \
   #-ERC GVCF
 
def haplotype_caller_region_bo(ref_genome_fasta, haplocaller_vcf, haplocaller_vcf_idx, bamout, individuals, chrom, start, end, includelist, completed):
  inputs = [ref_genome_fasta]
  outputs = [haplocaller_vcf, haplocaller_vcf_idx, bamout, completed]
  options = {"cores": 1 ,"memory": "16g","walltime":"160:00:00"}
  spec='''
  source activate prj1
  source /com/extra/java/8/load.sh
  source /com/extra/GATK/4/load.sh  

  java -Djava.io.tmpdir=/scratch/$PBS_JOBID/tmp \
     -Djava.awt.headless=true \
     -Xmx16g \
     -jar /com/extra/GATK/3.8/jar-bin/GenomeAnalysisTK.jar \
     -T HaplotypeCaller \
     -nct 1 \
     -R {ref_genome_fasta} \
     {includelist} \
     -A DepthPerSampleHC \
     -A Coverage \
     -A HaplotypeScore \
     -A StrandAlleleCountsBySample \
     -L {chrom}:{start}-{end} \
     -bamout {bamout} \
     -o {haplocaller_vcf}

  echo "Completed at "$(date) > {completed}'''.format(ref_genome_fasta=ref_genome_fasta, includelist = includelist, chrom = chrom, start=start, end=end, haplocaller_vcf=haplocaller_vcf, haplocaller_vcf_idx=haplocaller_vcf_idx, bamout=bamout, completed=completed)

  return inputs, outputs, options, spec

def combine_regions(split_file, species_haplocaller_vcf, specie, regions, refname, completed):
  inputs = [split_file]
  outputs = [species_haplocaller_vcf, completed]
  options = {"walltime":"10:00:00"}
  spec='''
  source activate prj1
  source /com/extra/vcftools/LATEST/load.sh

  vcf-concat $(cat {split_file} | fgrep -v chrX | fgrep -v chrUn | fgrep -v chrM | fgrep -v chrY | fgrep -v hap | sort -k1,1 -k2,2n | awk -v x=new_vcf_files_dir\/{specie}\/haplocaller.raw '{{OFS=".";print x,$1,$2,$3}}'|awk -v x=vcf '{{OFS="."; print $1,x}}') > {species_haplocaller_vcf}

  echo "Completed at "$(date) > {completed}'''.format(split_file=split_file, species_haplocaller_vcf=species_haplocaller_vcf, specie = specie, completed=completed)

  return inputs, outputs, options, spec

def get_readgroup(merged_bam_rg, read_group, completed):
  inputs = [merged_bam_rg]
  outputs = [read_group, completed]
  options = {"memory":"1g","walltime":"00:59:00"}
  spec='''
  source activate prj1

  cat {merged_bam_rg} | cut -f2 | cut -d: -f2 > {read_group}

  echo "Completed at "$(date) > {completed}'''.format(merged_bam_rg=merged_bam_rg, read_group=read_group, completed=completed)

  return inputs, outputs, options, spec

def split_bamout_region(read_group, bamout, bamout_individual_bam, bamout_individual_bai, completed):
  inputs = [read_group, bamout]
  outputs = [bamout_individual_bam, bamout_individual_bai, completed]
  options = {"memory":"4g","walltime":"10:00:00"}
  spec='''
  source activate prj1
  source /com/extra/samtools/LATEST/load.sh

  samtools view -R {read_group} -b -o {bamout_individual_bam} {bamout}
  samtools index {bamout_individual_bai}

  echo "Completed at "$(date) > {completed}'''.format(read_group = read_group, bamout = bamout, bamout_individual_bam = bamout_individual_bam, bamout_individual_bai = bamout_individual_bai, completed = completed)

  return inputs, outputs, options, spec

#inkluderer info om hvorvidt basen er repeat_masket

def get_family_coverage_region(filtered_bam_father, filtered_bam_mother, filtered_bam_child, family_coverage, completed):
  inputs = [filtered_bam_father, filtered_bam_mother, filtered_bam_child]
  outputs = [family_coverage, completed]
  options = {"memory":"1g","walltime":"10:00:00"}
  spec='''
  source activate prj1
  source /com/extra/samtools/LATEST/load.sh

  samtools depth {filtered_bam_father} {filtered_bam_mother} {filtered_bam_child} | awk '$3>=5 && $3<=150 && $4>=5 && $4<=150 && $5>=5 && $5<=150' | cut -f3- | /home/mtxellrb/MutationRate/scripts/table.py > {family_coverage}

  echo "Completed at "$(date) > {completed}'''.format(filtered_bam_father = filtered_bam_father, filtered_bam_mother = filtered_bam_mother, filtered_bam_child = filtered_bam_child, family_coverage = family_coverage, completed = completed)

  return inputs, outputs, options, spec

def combine_coverage_context_wr(input_files, new_family_coverage, chromosomes, child, minq, minQ, prev_completed, completed):
  inputs = [prev_completed]
  outputs = [new_family_coverage, completed]
  options = {"memory":"8g","walltime":"10:00:00"}
  spec='''
  source activate prj1

  input = $(echo {input_files}| tr -d "[]'"| tr "," " ")
  $input | cat | /home/mtxellrb/MutationRate/scripts/combine_tables.py > {new_family_coverage}

  echo "Completed at "$(date) > {completed}'''.format(input_files = input_files, new_family_coverage = new_family_coverage, completed = completed)

  return inputs, outputs, options, spec

def get_family_coverage_context_wr(filtered_bam_father, filtered_bai_father, filtered_bam_mother, filtered_bai_mother, filtered_bam_child, filtered_bai_child, family_coverage, minQ, minq, refname, chrom, completed):
  inputs = [filtered_bam_father, filtered_bai_father, filtered_bam_mother, filtered_bai_mother, filtered_bam_child, filtered_bai_child]
  outputs = [family_coverage, completed]
  options = {"walltime":"48:00:00"}
  spec='''
  source activate prj1
  source /com/extra/samtools/LATEST/load.sh
  source /com/extra/python/LATEST/load.sh

  samtools depth -Q{minQ} -q{minq} {filtered_bam_father} {filtered_bam_mother} {filtered_bam_child} -r{chrom}| awk '$3>=5 && $3<=200 && $4>=5 && $4<=200 && $5>=5 && $5<=200' | /home/mtxellrb/MutationRate/scripts/add_context_type_w_repeat.py {refname} | cut -f3- | /home/mtxellrb/MutationRate/scripts/table.py > {family_coverage}

  echo "Completed at "$(date) > {completed}'''.format(filtered_bam_father = filtered_bam_father, filter_bai_father = filtered_bam_father, filtered_bam_mother = filtered_bam_mother, filter_bai_mother = filtered_bam_mother, filtered_bam_child = filtered_bam_child, filter_bai_child = filtered_bam_child, family_coverage = family_coverage, minQ=minQ, minq=minq, refname = refname, chrom=chrom, completed = completed)

  return inputs, outputs, options, spec

######################################################################################
#C.PATHS AND VARIABLES
###################################################################################### 

pwd                 = "/home/mtxellrb/MutationRate/"
out                 = pwd+"out/"
compl_folder        = pwd+"out/COMPLETED/"

if not os.path.exists(out):
    os.makedirs(out)
if not os.path.exists(compl_folder):
    os.makedirs(compl_folder)

ref_genomes_dir     = pwd+"ref_genomes/"
fastq_files_dir     = pwd+"fastq_files/"

bam_files_dir       = pwd+"bam_files/"
merged_bams_dir_1   = pwd+"merged_bams/"
merged_bams_dir     = pwd+"actual_bam_files/"
filtered_bam_dir    = pwd+"filtered_bams/" 
new_gatk_files_dir  = pwd+"new_gatk_files/" 
split_files_dir     = pwd+"split_files/"
new_vcf_files_dir   = pwd+"new_vcf_files/"
read_group_dir      = pwd+"read_groups/"
bamout_files_dir    = pwd+"bamout_files/"
family_coverage_dir = pwd+"new_family_coverage_wr/"

if not os.path.exists(bam_files_dir):
    os.makedirs(bam_files_dir)
if not os.path.exists(merged_bams_dir_1):
    os.makedirs(merged_bams_dir_1)
if not os.path.exists(merged_bams_dir):
    os.makedirs(merged_bams_dir)
if not os.path.exists(filtered_bam_dir):
    os.makedirs(filtered_bam_dir)
if not os.path.exists(new_gatk_files_dir):
    os.makedirs(new_gatk_files_dir)
if not os.path.exists(split_files_dir):
    os.makedirs(split_files_dir)
if not os.path.exists(new_vcf_files_dir):
    os.makedirs(new_vcf_files_dir)
if not os.path.exists(read_group_dir):
    os.makedirs(read_group_dir)
if not os.path.exists(bamout_files_dir):
    os.makedirs(bamout_files_dir)
if not os.path.exists(family_coverage_dir):
    os.makedirs(family_coverage_dir)

CHIMP_REF     = 'panTro5' 
BONOBO_REF    = 'panPan2'
GORILLA_REF   = 'gorGor4'
ORANGUTAN_REF = 'ponAbe2'
GIBBON_REF    = 'nemLeu3'

BONOBOS = []
#CHIMPS     = ['Frits', 'Carolina', 'Carl', 'Simliki', 'ERR466113','ERR466114','ERR466115','ERR466116','ERR466117','ERR466118','ERR466119','ERR466120','ERR466121']
CHIMPS = ['Frits', 'Carolina', 'Carl']
#GORILLAS   = ['Banjo', 'Mimi', 'Mutasi', 'Mawenzi', 'Efata', 'Undi']
GORILLAS = []
#ORANGUTANS = ['Buschi', 'Moni', 'Masala', 'Farida', 'Ogan',  'Aris', 'Schubbi', 'Pongo']
ORANGUTANS = []
GIBBONS = []

CHIMP_FAMILIES = read_family_description('family_description/chimpanzees.txt')
BONOBO_FAMILIES = ""
GORILLA_FAMILIES = ""
ORANGUTAN_FAMILIES = ""
GIBBON_FAMILIES = ""
#BONOBO_FAMILIES = read_family_description('family_description/bonobos.txt')
#GORILLA_FAMILIES = read_family_description('family_description/gorillas.txt')
#ORANGUTAN_FAMILIES = read_family_description('family_description/orangutans.txt')
#GIBBON_FAMILIES = read_family_description('family_description/gibbons.txt')

#Order: (Child, Father, Mother)

ALL_APES = CHIMPS + BONOBOS + GORILLAS + ORANGUTANS + GIBBONS

#chimp_known = ['../faststorage/KnownVariants/pantro.snps.2012-11-25.autos.99vqsr.beagle_' + CHIMP_REF +'_lift.vcf',
#               '../faststorage/KnownVariants/gatk.allchimpanzee.2012-05-30.snps.chrX.vqsr99_' + CHIMP_REF + '_lift.vcf']

chimp_known = ["/home/mtxellrb/MutationRate/known_variants/Pan_troglodytes_RefPanTro5.vcf"]

#bonobo_known = [""]

#gorilla_known = ['../faststorage/KnownVariants/emitall-var.dropsamples.autos.vqsr99.2012-04-26.abfilter_' + GORILLA_REF  + '_lift_sorted.vcf']

#orang_known = ['../KnownVariants/abelii.2012-08-17.snps.dropsamples.autos.vqsr99.vcf',
#               '../KnownVariants/abelii.2012-08-17.snps.dropsamples.chrX.vcf',
#               '../KnownVariants/pygmaeus.2012-08-17.snps.dropsamples.autos.vqsr99.vcf',
#               '../KnownVariants/pygmaeus.2012-08-17.snps.dropsamples.chrX.vcf']

#gibbon_known = [""]

######################################################################################
#D.CODE
###################################################################################### 

#First we need to index the reference genomes with Picard
#refGenomes should be in ./ref-genomes

#REFS = [CHIMP_REF, BONOBO_REF, GORILLA_REF, ORANGUTAN_REF, GIBBON_REF]
REFS = [CHIMP_REF]

for ref in REFS:

  ref_genome_fasta = ref_genomes_dir+ref+".fa"
  ref_genome_fai = ref_genomes_dir+ref+".fa.fai"
  ref_genome_dict = ref_genomes_dir+ref+".dict"
  completed = compl_folder+"GetIndexRef_"+ref+".COMPLETED"

  gwf.target_from_template("GetIndexRef_"+ref,
          picard_index_reference(ref_genome_fasta = ref_genome_fasta,ref_genome_dict = ref_genome_dict,ref_genome_fai = ref_genome_fai, completed = completed))

  # Mapping fastq files to the reference

  ref_genome_am = ref_genomes_dir+ref+".am"
  ref_genome_ann = ref_genomes_dir+ref+".ann"
  ref_genome_pac = ref_genomes_dir+ref+".pac"
  ref_genome_sa = ref_genomes_dir+ref+".sa"
  completed = compl_folder+"GetIndexBwa_"+ref+".COMPLETED"

  gwf.target_from_template("GetIndexBwa_"+ref,
          get_index_bwa(ref = ref, ref_genome_fasta = ref_genome_fasta, ref_genome_am = ref_genome_am, ref_genome_ann = ref_genome_ann, ref_genome_pac = ref_genome_pac, ref_genome_sa = ref_genome_sa, completed = completed))

  reads = open(fastq_files_dir+'fastq_files_'+ref+'.txt').read().rstrip("\n").split("\n")
  
  name2fastq_pairs = {}
  for line in reads:
      name, fq1, fq2 = line.split(" ")
      if name.lower() == 'individual':
          continue
      if name not in name2fastq_pairs:
          name2fastq_pairs[name] = []
      name2fastq_pairs[name].append((fq1, fq2))
    
  for name in name2fastq_pairs:
    bam_files = []
    for i in range(len(name2fastq_pairs[name])):
      read_1, read_2 = name2fastq_pairs[name][i]
      bam = bam_files_dir+name+"_"+str(i)+".bam"
      bam_files.append(bam)

      prev_completed = compl_folder+"GetIndexBwa_"+ref+".COMPLETED"
      completed = compl_folder+"Mapping_"+name+"_"+str(i)+"_"+ref+".COMPLETED"
    
      gwf.target_from_template("Mapping_"+name+"_"+str(i)+"_"+ref,
              map_bwa(ref = ref, read_1 = read_1, read_2 = read_2, ref_genome_am = ref_genome_am, ref_genome_ann = ref_genome_ann, ref_genome_pac = ref_genome_pac, ref_genome_sa = ref_genome_sa, bam = bam, prev_completed = prev_completed, completed = completed))

    prev_completed = compl_folder+"Mapping_"+name+"_"+str(len(name2fastq_pairs[name])-1)+"_"+ref+".COMPLETED"
    inputbams = " ".join(bam_files)
    merged_bam = merged_bams_dir_1+name+".bam"
    merged_bam_rg = merged_bams_dir_1+name+".rg.txt"
    completed = compl_folder+"MergeBams_"+name+"_"+ref+".COMPLETED"    

    gwf.target_from_template("MergeBams_"+name+"_"+ref,
            merge_bams(inputbams = inputbams, merged_bam = merged_bam, merged_bam_rg = merged_bam_rg, name = name, prev_completed = prev_completed, completed = completed))

# Use GATK to process bam files before calling
#To work with the bam files in GATK we first need them indexed...

for individual in ALL_APES:

  merged_bam = merged_bams_dir+individual+".bam"
  merged_bai = merged_bams_dir+individual+".bam.bai"
  completed = compl_folder+"GetIndexMergedBam_"+individual+".COMPLETED"

  gwf.target_from_template("GetIndexMergedBam_"+individual,
          samtools_index_bam(merged_bam = merged_bam, merged_bai = merged_bai, completed = completed))

  filtered_bam = filtered_bam_dir+individual+".bam"
  filtered_bai = filtered_bam_dir+individual+".bam.bai"
  completed = compl_folder+"GetIndexFilteredBam_"+individual+".COMPLETED"

  gwf.target_from_template("GetIndexFilteredBam_"+individual,
          samtools_index_bam(filtered_bam = filtered_bam, filtered_bai = filtered_bai, completed = completed))

#Next, identify regions that needs to be locally re-aligned.
# intermediary output files from gatk will be put in ./gatk_files/

#Realign the bam files:

for individual in CHIMPS:

  ref_genome_fasta = ref_genomes_dir+CHIMP_REF+".fa"
  ref_genome_dict = ref_genomes_dir+CHIMP_REF+".dict"
  merged_bam = merged_bams_dir+individual+".bam"
  merged_bai = merged_bams_dir+individual+".bam.bai"
  recalibration_report = new_gatk_files_dir+individual+".recalibration_report.grp"
  completed = compl_folder+"CalcRecInfo_"+individual+".COMPLETED"
  
  prev_completed = compl_folder+"GetIndexMergedBam_"+individual+".COMPLETED"

  knownlist = ' '.join(x for x in chimp_known)

  gwf.target_from_template("CalcRecInfo_"+individual,
        calc_recalibrate_info(ref_genome_fasta=ref_genome_fasta, ref_genome_dict=ref_genome_dict, merged_bam=merged_bam, merged_bai=merged_bai, knownlist=knownlist, recalibration_report=recalibration_report, prev_completed = prev_completed, completed=completed))

  recalibrated_bam = new_gatk_files_dir+individual+".recalibrated.bam"
  recalibrated_bai = new_gatk_files_dir+individual+".recalibrated.bam.bai"
  completed = compl_folder+"GetRecalBam_"+individual+".COMPLETED"

  prev_completed = compl_folder+"GetIndexMergedBam_"+individual+".COMPLETED"

  gwf.target_from_template("GetRecalBam_"+individual,
        recalibrate(ref_genome_fasta = ref_genome_fasta, ref_genome_dict = ref_genome_dict, merged_bam = merged_bam, merged_bai = merged_bai, recalibration_report = recalibration_report, recalibrated_bam = recalibrated_bam, recalibrated_bai = recalibrated_bai, prev_completed = prev_completed, completed = completed))

  filtered_bam = filtered_bam_dir+individual+".bam"
  completed = compl_folder+"FilterBam_"+individual+".COMPLETED"

  prev_completed = compl_folder+"GetRecalBam_"+individual+".COMPLETED"

  gwf.target_from_template("FilterBam_"+individual,
        filter_bam_files(recalibrated_bam = recalibrated_bam, recalibrated_bai = recalibrated_bai, ref_genome_fasta = ref_genome_fasta, filtered_bam = filtered_bam, prev_completed = prev_completed, completed = completed))

for individual in BONOBOS:

  ref_genome_fasta = ref_genomes_dir+BONOBO_REF+".fa"
  ref_genome_dict = ref_genomes_dir+BONOBO_REF+".dict"
  merged_bam = merged_bams_dir+individual+".bam"
  merged_bai = merged_bams_dir+individual+".bam.bai"
  recalibration_report = new_gatk_files_dir+individual+".recalibration_report.grp"
  completed = compl_folder+"CalcRecInfo_"+individual+".COMPLETED"

  prev_completed = compl_folder+"GetIndexMergedBam_"+individual+".COMPLETED"

  knownlist = ' '.join(x for x in bonobo_known)

  gwf.target_from_template("CalcRecInfo_"+individual,
        calc_recalibrate_info(ref_genome_fasta=ref_genome_fasta, ref_genome_dict=ref_genome_dict, merged_bam=merged_bam, merged_bai=merged_bai, knownlist=chimp_known, recalibration_report=recalibration_report, prev_completed=prev_completed, completed=completed))

  recalibrated_bam = new_gatk_files_dir+individual+".recalibrated.bam"
  recalibrated_bai = new_gatk_files_dir+individual+".recalibrated.bam.bai"
  completed = compl_folder+"GetRecalBam_"+individual+".COMPLETED"

  prev_completed = compl_folder+"GetIndexMergedBam_"+individual+".COMPLETED"

  gwf.target_from_template("GetRecalBam_"+individual,
        recalibrate(ref_genome_fasta = ref_genome_fasta, ref_genome_dict = ref_genome_dict, merged_bam = merged_bam, merged_bai = merged_bai, recalibration_report = recalibration_report, recalibrated_bam = recalibrated_bam, recalibrated_bai = recalibrated_bai, prev_completed = prev_completed, completed = completed))

  filtered_bam = filtered_bam_dir+individual+".bam"
  completed = compl_folder+"FilterBam_"+individual+".COMPLETED"

  prev_completed = compl_folder+"GetRecalBam_"+individual+".COMPLETED"

  gwf.target_from_template("FilterBam_"+individual,
        filter_bam_files(recalibrated_bam = recalibrated_bam, recalibrated_bai = recalibrated_bai, ref_genome_fasta = ref_genome_fasta, filtered_bam = filtered_bam, prev_completed = prev_completed, completed = completed))

for individual in GORILLAS:

  ref_genome_fasta = ref_genomes_dir+GORILLA_REF+".fa"
  ref_genome_dict = ref_genomes_dir+GORILLA_REF+".dict"
  merged_bam = merged_bams_dir+individual+".bam"
  merged_bai = merged_bams_dir+individual+".bam.bai"
  recalibration_report = new_gatk_files_dir+individual+".recalibration_report.grp"
  completed = compl_folder+"CalcRecInfo_"+individual+".COMPLETED"

  prev_completed = compl_folder+"GetIndexMergedBam_"+individual+".COMPLETED"

  knownlist = ' '.join(x for x in gorilla_known)

  gwf.target_from_template("CalcRecInfo_"+individual,
        calc_recalibrate_info(ref_genome_fasta=ref_genome_fasta, ref_genome_dict=ref_genome_dict, merged_bam=merged_bam, merged_bai=merged_bai, knownlist=gorilla_known, recalibration_report=recalibration_report, prev_completed = prev_completed, completed=completed))

  recalibrated_bam = new_gatk_files_dir+individual+".recalibrated.bam"
  recalibrated_bai = new_gatk_files_dir+individual+".recalibrated.bam.bai"
  completed = compl_folder+"GetRecalBam_"+individual+".COMPLETED"

  prev_completed = compl_folder+"GetIndexMergedBam_"+individual+".COMPLETED"

  gwf.target_from_template("GetRecalBam_"+individual,
        recalibrate(ref_genome_fasta = ref_genome_fasta, ref_genome_dict = ref_genome_dict, merged_bam = merged_bam, merged_bai = merged_bai, recalibration_report = recalibration_report, recalibrated_bam = recalibrated_bam, recalibrated_bai = recalibrated_bai, prev_completed = prev_completed, completed = completed))

  filtered_bam = filtered_bam_dir+individual+".bam"
  completed = compl_folder+"FilterBam_"+individual+".COMPLETED"

  prev_completed = compl_folder+"GetRecalBam_"+individual+".COMPLETED"

  gwf.target_from_template("FilterBam_"+individual,
        filter_bam_files(recalibrated_bam = recalibrated_bam, recalibrated_bai = recalibrated_bai, ref_genome_fasta = ref_genome_fasta, filtered_bam = filtered_bam, prev_completed = prev_completed, completed = completed))

for individual in ORANGUTANS:

  ref_genome_fasta = ref_genomes_dir+ORANGUTAN_REF+".fa"
  ref_genome_dict = ref_genomes_dir+ORANGUTAN_REF+".dict"
  merged_bam = merged_bams_dir+individual+".bam"
  merged_bai = merged_bams_dir+individual+".bam.bai"
  recalibration_report = new_gatk_files_dir+individual+".recalibration_report.grp"
  completed = compl_folder+"CalcRecInfo_"+individual+".COMPLETED"

  prev_completed = compl_folder+"GetIndexMergedBam_"+individual+".COMPLETED"

  knownlist = ' '.join(x for x in orangutan_known)

  gwf.target_from_template("CalcRecInfo_"+individual,
        calc_recalibrate_info(ref_genome_fasta=ref_genome_fasta, ref_genome_dict=ref_genome_dict, merged_bam=merged_bam, merged_bai=merged_bai, knownlist=orang_known, recalibration_report=recalibration_report, prev_completed = prev_completed, completed=completed))

  recalibrated_bam = new_gatk_files_dir+individual+".recalibrated.bam"
  recalibrated_bai = new_gatk_files_dir+individual+".recalibrated.bam.bai"
  completed = compl_folder+"GetRecalBam_"+individual+".COMPLETED"

  prev_completed = compl_folder+"GetIndexMergedBam_"+individual+".COMPLETED"

  gwf.target_from_template("GetRecalBam_"+individual,
        recalibrate(ref_genome_fasta = ref_genome_fasta, ref_genome_dict = ref_genomes_dict, merged_bam = merged_bam, merged_bai = merged_bai, recalibration_report = recalibration_report, recalibrated_bam = recalibrated_bam, recalibrated_bai = recalibrated_bai, prev_completed = prev_completed, completed = completed))

  filtered_bam = filtered_bam_dir+individual+".bam"
  completed = compl_folder+"FilterBam_"+individual+".COMPLETED"

  prev_completed = compl_folder+"GetRecalBam_"+individual+".COMPLETED"

  gwf.target_from_template("FilterBam_"+individual,
        filter_bam_files(recalibrated_bam = recalibrated_bam, recalibrated_bai = recalibrated_bai, ref_genome_fasta = ref_genome_fasta, filtered_bam = filtered_bam, prev_completed = prev_completed, completed = completed))

for individual in GIBBONS:

  ref_genome_fasta = ref_genomes_dir+GIBBON_REF+".fa"
  ref_genome_dict = ref_genomes_dir+GIBBON_REF+".dict"
  merged_bam = merged_bams_dir+individual+".bam"
  merged_bai = merged_bams_dir+individual+".bam.bai"
  recalibration_report = new_gatk_files_dir+individual+".recalibration_report.grp"
  completed = compl_folder+"CalcRecInfo_"+individual+".COMPLETED"

  prev_completed = compl_folder+"GetIndexMergedBam_"+individual+".COMPLETED"

  knownlist = ' '.join(x for x in gibbon_known)

  gwf.target_from_template("CalcRecInfo_"+individual,
        calc_recalibrate_info(ref_genome_fasta=ref_genome_fasta, ref_genome_dict=ref_genome_dict, merged_bam=merged_bam, merged_bai=merged_bai, knownlist=chimp_known, recalibration_report=recalibration_report, prev_completed = prev_completed, completed=completed))

  recalibrated_bam = new_gatk_files_dir+individual+".recalibrated.bam"
  recalibrated_bai = new_gatk_files_dir+individual+".recalibrated.bam.bai"
  completed = compl_folder+"GetRecalBam_"+individual+".COMPLETED"

  prev_completed = compl_folder+"GetIndexMergedBam_"+individual+".COMPLETED"

  gwf.target_from_template("GetRecalBam_"+individual,
        recalibrate(ref_genome_fasta = ref_genome_fasta, ref_genome_dict = ref_genomes_dict, merged_bam = merged_bam, merged_bai = merged_bai, recalibration_report = recalibration_report, recalibrated_bam = recalibrated_bam, recalibrated_bai = recalibrated_bai, prev_completed = prev_completed, completed = completed))

  filtered_bam = filtered_bam_dir+individual+".bam"
  completed = compl_folder+"FilterBam_"+individual+".COMPLETED"

  prev_completed = compl_folder+"GetRecalBam_"+individual+".COMPLETED"

  gwf.target_from_template("FilterBam_"+individual,
        filter_bam_files(recalibrated_bam = recalibrated_bam, recalibrated_bai = recalibrated_bai, ref_genome_fasta = ref_genome_fasta, filtered_bam = filtered_bam, prev_completed = prev_completed, completed = completed))

#Call variants using haplotype caller

def get_chromosomes(refname):
    f = open(ref_genomes_dir+ refname + '.fa.fai')
    chroms = []
    for line in f:
        chrom = line.split()[0]
        #if 'random' not in chrom and
        if 'Un' not in chrom:
            chroms.append(chrom)
    f.close()
    return chroms

def is_autosome(chrom):
    return (not chrom.startswith('chrX') and
            not chrom.startswith('chrY') and
            not chrom.startswith('chrUn') and
            not chrom.startswith('chrM'))

# Regions used to split: Find that file

def get_regions(refname):
    f = open(split_files_dir + refname + '_split_1000_1000000.txt')
    L = []
    for line in f:
        chrom, start, end = line.split()
        L.append((chrom, int(start), int(end)))
    f.close()
    return L

CHIMP_CHROMOSOMES = get_chromosomes(CHIMP_REF)
BONOBO_CHROMOSOMES = []
GORILLA_CHROMOSOMES = []
ORANGUTAN_CHROMOSOMES = []
GIBBON_CHROMOSOMES = []
#BONOBO_CHROMOSOMES = get_chromosomes(BONOBO_REF)
#GORILLA_CHROMOSOMES = get_chromosomes(GORILLA_REF)
#ORANGUTAN_CHROMOSOMES = get_chromosomes(ORANGUTAN_REF)
#GIBBON_CHROMOSOMES = get_chromosomes(GIBBON_REF)

CHIMP_REGIONS = get_regions(CHIMP_REF)
BONOBO_REGIONS = []
GORILLA_REGIONS = []
ORANGUTAN_REGIONS = []
GIBBON_REGIONS = []
#BONOBO_REGIONS = get_regions(BONOBO_REF)
#GORILLA_REGIONS = get_regions(GORILLA_REF)
#ORANGUTAN_REGIONS = get_regions(ORANGUTAN_REF)
#GIBBON_REGIONS = get_regions(GIBBON_REF)

CHIMP_AUTOSOME_REGIONS = [(chrom, start, end) for (chrom, start, end) in CHIMP_REGIONS if not (chrom.startswith('chrX') or chrom.startswith('chrY') or chrom.startswith('chrM'))]
CHIMP_X_REGIONS = [(chrom, start, end) for (chrom, start, end) in CHIMP_REGIONS if chrom.startswith('chrX')]

BONOBO_AUTOSOME_REGIONS = []
BONOBO_X_REGIONS = []
#BONOBO_AUTOSOME_REGIONS = [(chrom, start, end) for (chrom, start, end) in BONOBO_REGIONS if not (chrom.startswith('chrX') or chrom.startswith('chrY') or chrom.startswith('chrM'))]
#BONOBO_X_REGIONS = [(chrom, start, end) for (chrom, start, end) in BONOBO_REGIONS if chrom.startswith('chrX')]

GORILLA_AUTOSOME_REGIONS = []
GORILLA_X_REGIONS = []
#GORILLA_AUTOSOME_REGIONS = [(chrom, start, end) for (chrom, start, end) in GORILLA_REGIONS if not (chrom.startswith('chrX') or chrom.startswith('chrY') or chrom.startswith('chrM'))]
#GORILLA_X_REGIONS = [(chrom, start, end) for (chrom, start, end) in GORILLA_REGIONS if chrom.startswith('chrX')]

ORANGUTAN_AUTOSOME_REGIONS = []
ORANGUTAN_X_REGIONS = []
#ORANGUTAN_AUTOSOME_REGIONS = [(chrom, start, end) for (chrom, start, end) in ORANGUTAN_REGIONS if not (chrom.startswith('chrX') or chrom.startswith('chrY') or chrom.startswith('chrM') or 'hap' in chrom)]
#ORANGUTAN_X_REGIONS = [(chrom, start, end) for (chrom, start, end) in ORANGUTAN_REGIONS if chrom.startswith('chrX')]

GIBBON_AUTOSOME_REGIONS = []
GIBBON_X_REGIONS = []

CHIMP_AUTOSOMES = [x for x in CHIMP_CHROMOSOMES if is_autosome(x)]
BONOBO_AUTOSOMES = []
GORILLA_AUTOSOMES = []
ORANGUTAN_AUTOSOMES = []
GIBBON_AUTOSOMES = []
#BONOBO_AUTOSOMES = [x for x in BONOBO_CHROMOSOMES if is_autosome(x)]
#GORILLA_AUTOSOMES = [x for x in GORILLA_CHROMOSOMES if is_autosome(x)]
#ORANGUTAN_AUTOSOMES = [x for x in ORANGUTAN_CHROMOSOMES if is_autosome(x)]
#GIBBON_AUTOSOMES = [x for x in GIBBON_CHROMOSOMES if is_autosome(x)]

for chrom, start, end in CHIMP_REGIONS:

  chrom = str(chrom)
  start = str(start)
  end = str(end)

  individuals = CHIMPS
  ref_genome_fasta = ref_genomes_dir+CHIMP_REF+".fa"
  specie = "chimpanzees"

  haplocaller_out_dir = new_vcf_files_dir+specie+"/"
  bamout_dir = new_gatk_files_dir+"bamout/"+specie+"/"

  if not os.path.exists(haplocaller_out_dir):
    os.makedirs(haplocaller_out_dir)

  if not os.path.exists(bamout_dir):
    os.makedirs(bamout_dir)

  haplocaller_vcf = haplocaller_out_dir+"haplocaller.raw."+chrom+"."+start+"."+end+".vcf"
  haplocaller_vcf_idx = haplocaller_out_dir+"haplocaller.raw."+chrom+"."+start+"."+end+".vcf.idx"
  bamout = bamout_dir+chrom+"."+start+"."+end+".bam"
  completed = compl_folder+"HaplotypeCaller_"+chrom+"_"+start+"_"+end+".COMPLETED"

  prev_completed = compl_folder+"GetRecalBam_"+individual+".COMPLETED"

  bamfiles = [new_gatk_files_dir+individual+'.recalibrated.bam' for individual in individuals]
  includelist = ' '.join('-I ' + x for x in bamfiles)

  gwf.target_from_template("HaplotypeCaller_"+chrom+"_"+start+"_"+end,
          haplotype_caller_region_bo(ref_genome_fasta=ref_genome_fasta, haplocaller_vcf=haplocaller_vcf, haplocaller_vcf_idx=haplocaller_vcf_idx, bamout=bamout, individuals=individuals, chrom=chrom, start=start, end=end, includelist=includelist, prev_completed = prev_completed, completed=completed))

  split_file = split_files_dir + CHIMP_REF + "_split_1000_1000000.txt"
  species_haplocaller_vcf = new_vcf_files_dir+specie+"_"+chrom+"_"+start+"_"+end+".haplocaller.raw.auto.vcf"
  completed = compl_folder+"CombineRegions_Chimps_"+chrom+"_"+start+"_"+end+".COMPLETED"

  prev_completed = compl_folder+"HaplotypeCaller_"+chrom+"_"+start+"_"+end+".COMPLETED"

  #vcf_files = "".join([new_vcf_files_dir + specie + '/haplocaller.raw.' + chrom + '.' + str(start) + '.' + str(end) + '.vcf' for (chrom,start,end) in CHIMP_REGIONS])

  gwf.target_from_template("CombineRegions_Chimps_"+chrom+"_"+start+"_"+end,
          combine_regions(refname=CHIMP_REF, specie=specie, regions=CHIMP_REGIONS, split_file = split_file, species_haplocaller_vcf = species_haplocaller_vcf, prev_completed = prev_completed, completed = completed))

for chrom, start, end in BONOBO_REGIONS:

  chrom = str(chrom)
  start = str(start)
  end = str(end)

  individuals = BONOBOS
  ref_genome_fasta = ref_genomes_dir+BONOBO_REF+".fa"
  specie = "bonobos"

  haplocaller_out_dir = new_vcf_files_dir+specie+"/"
  bamout_dir = new_gatk_files_dir+"bamout/"+specie+"/"

  if not os.path.exists(haplocaller_out_dir):
    os.makedirs(haplocaller_out_dir)

  if not os.path.exists(bamout_dir):
    os.makedirs(bamout_dir)

  haplocaller_vcf = haplocaller_out_dir+"haplocaller.raw."+chrom+"."+start+"."+end+".vcf"
  haplocaller_vcf_idx = haplocaller_out_dir+"haplocaller.raw."+chrom+"."+start+"."+end+".vcf.idx"
  bamout = bamout_dir+chrom+"."+start+"."+end+".bam"
  completed = compl_folder+"HaplotypeCaller_"+chrom+"_"+start+"_"+end+".COMPLETED"

  prev_completed = compl_folder+"GetRecalBam_"+individual+".COMPLETED"

  bamfiles = [new_gatk_files_dir+individual+'.recalibrated.bam' for individual in individuals]
  includelist = ' '.join('-I ' + x for x in bamfiles)

  gwf.target_from_template("HaplotypeCaller_"+chrom+"_"+start+"_"+end,
          haplotype_caller_region_bo(ref_genome_fasta=ref_genome_fasta, haplocaller_vcf=haplocaller_vcf, haplocaller_vcf_idx=haplocaller_vcf_idx, bamout=bamout, individuals=individuals, chrom=chrom, start=start, end=end, includelist=includelist, prev_completed = prev_completed, completed=completed))

  split_file = split_files_dir + BONOBO_REF + "_split_1000_1000000.txt"
  species_haplocaller_vcf = new_vcf_files_dir+specie+"_"+chrom+"_"+start+"_"+end+".haplocaller.raw.auto.vcf"
  completed = compl_folder+"CombineRegions_Bonobos_"+chrom+"_"+start+"_"+end+".COMPLETED"

  prev_completed = compl_folder+"HaplotypeCaller_"+chrom+"_"+start+"_"+end+".COMPLETED"

  gwf.target_from_template("CombineRegions_Bonobos_"+chrom+"_"+start+"_"+end,
          combine_regions(refname=BONOBO_REF, specie=specie, regions=BONOBO_REGIONS, split_file = split_file, species_haplocaller_vcf = species_haplocaller_vcf, prev_completed = prev_completed, completed = completed))

for chrom, start, end in GORILLA_REGIONS:

  chrom = str(chrom)
  start = str(start)
  end = str(end)

  individuals = GORILLAS
  ref_genome_fasta = ref_genomes_dir+GORILLA_REF+".fa"
  specie = "gorillas"

  haplocaller_out_dir = new_vcf_files_dir+specie+"/"
  bamout_dir = new_gatk_files_dir+"bamout/"+specie+"/"

  if not os.path.exists(haplocaller_out_dir):
    os.makedirs(haplocaller_out_dir)

  if not os.path.exists(bamout_dir):
    os.makedirs(bamout_dir)

  haplocaller_vcf = haplocaller_out_dir+"haplocaller.raw."+chrom+"."+start+"."+end+".vcf"
  haplocaller_vcf_idx = haplocaller_out_dir+"haplocaller.raw."+chrom+"."+start+"."+end+".vcf.idx"
  bamout = bamout_dir+chrom+"."+start+"."+end+".bam"
  completed = compl_folder+"HaplotypeCaller_"+chrom+"_"+start+"_"+end+".COMPLETED"

  prev_completed = compl_folder+"GetRecalBam_"+individual+".COMPLETED"

  bamfiles = [new_gatk_files_dir+individual+'.recalibrated.bam' for individual in individuals]
  includelist = ' '.join('-I ' + x for x in bamfiles)

  gwf.target_from_template("HaplotypeCaller_"+chrom+"_"+start+"_"+end,
          haplotype_caller_region_bo(ref_genome_fasta=ref_genome_fasta, haplocaller_vcf=haplocaller_vcf, haplocaller_vcf_idx=haplocaller_vcf_idx, bamout=bamout, individuals=individuals, chrom=chrom, start=start, end=end, includelist=includelist, prev_completed = prev_completed, completed=completed))

  split_file = split_files_dir + GORILLA_REF + "_split_1000_1000000.txt"
  species_haplocaller_vcf = new_vcf_files_dir+specie+"_"+chrom+"_"+start+"_"+end+".haplocaller.raw.auto.vcf"
  completed = compl_folder+"CombineRegions_Gorillas_"+chrom+"_"+start+"_"+end+".COMPLETED"

  prev_completed = compl_folder+"HaplotypeCaller_"+chrom+"_"+start+"_"+end+".COMPLETED"

  gwf.target_from_template("CombineRegions_Gorillas_"+chrom+"_"+start+"_"+end,
          combine_regions(refname=GORILLA_REF, specie=specie, regions=GORILLAS_REGIONS, split_file = split_file, species_haplocaller_vcf = species_haplocaller_vcf, prev_completed = prev_completed, completed = completed))

for chrom, start, end in ORANGUTAN_AUTOSOME_REGIONS:

  chrom = str(chrom)
  start = str(start)
  end = str(end)

  individuals = ORANGUTANS
  ref_genome_fasta = ref_genomes_dir+ORANGUTAN_REF+".fa"
  specie = "orangutans"

  haplocaller_out_dir = new_vcf_files_dir+specie+"/"
  bamout_dir = new_gatk_files_dir+"bamout/"+specie+"/"

  if not os.path.exists(haplocaller_out_dir):
    os.makedirs(haplocaller_out_dir)

  if not os.path.exists(bamout_dir):
    os.makedirs(bamout_dir)

  haplocaller_vcf = haplocaller_out_dir+"haplocaller.raw."+chrom+"."+start+"."+end+".vcf"
  haplocaller_vcf_idx = haplocaller_out_dir+"haplocaller.raw."+chrom+"."+start+"."+end+".vcf.idx"
  bamout = bamout_dir+chrom+"."+start+"."+end+".bam"
  completed = compl_folder+"HaplotypeCaller_"+chrom+"_"+start+"_"+end+".COMPLETED"
  
  prev_completed = compl_folder+"GetRecalBam_"+individual+".COMPLETED"

  bamfiles = [new_gatk_files_dir+individual+'.recalibrated.bam' for individual in individuals]
  includelist = ' '.join('-I ' + x for x in bamfiles)

  gwf.target_from_template("HaplotypeCaller_"+chrom+"_"+start+"_"+end,
          haplotype_caller_region_bo(ref_genome_fasta=ref_genome_fasta, haplocaller_vcf=haplocaller_vcf, haplocaller_vcf_idx=haplocaller_vcf_idx, bamout=bamout, individuals=individuals, chrom=chrom, start=start, end=end, includelist=includelist, prev_completed = prev_completed, completed=completed))

  split_file = split_files_dir + ORANGUTAN_REF + "_split_1000_1000000.txt"
  species_haplocaller_vcf = new_vcf_files_dir+specie+"_"+chrom+"_"+start+"_"+end+".haplocaller.raw.auto.vcf"
  completed = compl_folder+"CombineRegions_Orangutans_"+chrom+"_"+start+"_"+end+"_"+".COMPLETED"

  prev_completed = compl_folder+"HaplotypeCaller_"+chrom+"_"+start+"_"+end+".COMPLETED"

  gwf.target_from_template("CombineRegions_Orangutans_"+chrom+"_"+start+"_"+end,
          combine_regions(refname=ORANGUTAN_REF, specie=specie, regions=ORANGUTAN_REGIONS, split_file = split_file, species_haplocaller_vcf = species_haplocaller_vcf, prev_completed = prev_completed, completed = completed))

for chrom, start, end in GIBBON_AUTOSOME_REGIONS:

  chrom = str(chrom)
  start = str(start)
  end = str(end)

  individuals = GIBBONS
  ref_genome_fasta = ref_genomes_dir+GIBBON_REF+".fa"
  specie = "gibbons"

  haplocaller_out_dir = new_vcf_files_dir+specie+"/"
  bamout_dir = new_gatk_files_dir+"bamout/"+specie+"/"

  if not os.path.exists(haplocaller_out_dir):
    os.makedirs(haplocaller_out_dir)

  if not os.path.exists(bamout_dir):
    os.makedirs(bamout_dir)

  haplocaller_vcf = haplocaller_out_dir+"haplocaller.raw."+chrom+"."+start+"."+end+".vcf"
  haplocaller_vcf_idx = haplocaller_out_dir+"haplocaller.raw."+chrom+"."+start+"."+end+".vcf.idx"
  bamout = bamout_dir+chrom+"."+start+"."+end+".bam"
  completed = compl_folder+"HaplotypeCaller_"+chrom+"_"+start+"_"+end+".COMPLETED"

  prev_completed = compl_folder+"GetRecalBam_"+individual+".COMPLETED"

  bamfiles = [new_gatk_files_dir+individual+'.recalibrated.bam' for individual in individuals]
  includelist = ' '.join('-I ' + x for x in bamfiles)

  gwf.target_from_template("HaplotypeCaller_"+chrom+"_"+start+"_"+end,
          haplotype_caller_region_bo(ref_genome_fasta=ref_genome_fasta, haplocaller_vcf=haplocaller_vcf, haplocaller_vcf_idx=haplocaller_vcf_idx, bamout=bamout, individuals=individuals, chrom=chrom, start=start, end=end, includelist=includelist, prev_completed = prev_completed, completed=completed))

  split_file = split_files_dir + GIBBON_REF + "_split_1000_1000000.txt"
  species_haplocaller_vcf = new_vcf_files_dir+specie+"_"+chrom+"_"+start+"_"+end+".haplocaller.raw.auto.vcf"
  completed = compl_folder+"CombineRegions_Gibbons_"+chrom+"_"+start+"_"+end+"_"+".COMPLETED"

  prev_completed = compl_folder+"HaplotypeCaller_"+chrom+"_"+start+"_"+end+".COMPLETED"

  gwf.target_from_template("CombineRegions_Gibbons_"+chrom+"_"+start+"_"+end,
          combine_regions(refname=GIBBON_REF, specie=specie, regions=GIBBON_REGIONS, split_file = split_file, species_haplocaller_vcf = species_haplocaller_vcf, prev_completed = prev_completed, completed = completed))

for individual in ALL_APES:

  merged_bam_rg = merged_bams_dir+individual+".rg.txt"
  read_group = read_group_dir+individual+".rg.txt"
  completed = compl_folder+"GetReadGroup_"+individual+".COMPLETED"

  prev_completed = compl_folder+"GetIndexMergedBam_"+individual+".COMPLETED"

  gwf.target_from_template("GetReadGroup_"+individual,
          get_readgroup(merged_bam_rg = merged_bam_rg, read_group = read_group, prev_completed = prev_completed, completed = completed))

for individual in set(itertools.chain(*CHIMP_FAMILIES)):

  for chrom, start, end in CHIMP_REGIONS:

    chrom = str(chrom)
    start = str(start)
    end = str(end)

    specie = "chimpanzees"
    bamout_dir = new_gatk_files_dir+"bamout/"+specie+"/"
    bamout = bamout_dir+chrom+"."+start+"."+end+".bam"
    read_group = read_group_dir+individual+".rg.txt"

    bamout_individual_dir = bamout_files_dir+individual+"/"

    if not os.path.exists(bamout_individual_dir):
      os.makedirs(bamout_individual_dir)

    bamout_individual_bam = bamout_files_dir+individual+"/"+chrom+"."+start+"."+end+".bam"
    bamout_individual_bai = bamout_files_dir+individual+"/"+chrom+"."+start+"."+end+".bai"
    completed = compl_folder+"SplitBamOutRegion_"+individual+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

    prev_completed = compl_folder+"GetReadGroup_"+individual+".COMPLETED"

    gwf.target_from_template("SplitBamOutRegion_"+individual+"_"+chrom+"_"+start+"_"+end,
            split_bamout_region(read_group = read_group, bamout_individual_bam = bamout_individual_bam, bamout_individual_bai = bamout_individual_bai, bamout = bamout, prev_completed = prev_completed, completed = completed))

for individual in set(itertools.chain(*BONOBO_FAMILIES)):
  for chrom, start, end in BONOBO_REGIONS:

    chrom = str(chrom)
    start = str(start)
    end = str(end)

    specie = "bonobos"
    bamout_dir = new_gatk_files_dir+"bamout/"+specie+"/"
    bamout = bamout_dir+chrom+"."+start+"."+end+".bam"
    read_group = read_group_dir+individual+".rg.txt"

    bamout_individual_dir = bamout_files_dir+individual+"/"

    if not os.path.exists(bamout_individual_dir):
      os.makedirs(bamout_individual_dir)

    bamout_individual_bam = bamout_individual_dir+chrom+"."+start+"."+end+".bam"
    bamout_individual_bai = bamout_individual_dir+chrom+"."+start+"."+end+".bai"
    completed = compl_folder+"SplitBamOutRegion_"+individual+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

    prev_completed = compl_folder+"GetReadGroup_"+individual+".COMPLETED"

    gwf.target_from_template("SplitBamOutRegion_"+individual+"_"+chrom+"_"+start+"_"+end,
            split_bamout_region(read_group = read_group, bamout_individual_bam = bamout_individual_bam, bamout_individual_bai = bamout_individual_bai, bamout = bamout, prev_completed = prev_completed, completed = completed))

for individual in set(itertools.chain(*GORILLA_FAMILIES)):
  for chrom, start, end in GORILLA_REGIONS:

    chrom = str(chrom)
    start = str(start)
    end = str(end)

    specie = "gorillas"
    bamout_dir = new_gatk_files_dir+"bamout/"+specie+"/"
    bamout = bamout_dir+chrom+"."+start+"."+end+".bam"
    read_group = read_group_dir+individual+".rg.txt"

    bamout_individual_dir = bamout_files_dir+individual+"/"

    if not os.path.exists(bamout_individual_dir):
      os.makedirs(bamout_individual_dir)

    bamout_individual_bam = bamout_individual_dir+chrom+"."+start+"."+end+".bam"
    bamout_individual_bai = bamout_individual_dir+chrom+"."+start+"."+end+".bai"
    completed = compl_folder+"SplitBamOutRegion_"+individual+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

    prev_completed = compl_folder+"GetReadGroup_"+individual+".COMPLETED"

    gwf.target_from_template("SplitBamOutRegion_"+individual+"_"+chrom+"_"+start+"_"+end,
            split_bamout_region(read_group = read_group, bamout_individual_bam = bamout_individual_bam, bamout_individual_bai = bamout_individual_bai,bamout=bamout, prev_completed = prev_completed, completed = completed))

for individual in set(itertools.chain(*ORANGUTAN_FAMILIES)):
  for chrom, start, end in ORANGUTAN_AUTOSOME_REGIONS:

    chrom = str(chrom)
    start = str(start)
    end = str(end)

    specie = "orangutans"
    bamout_dir = new_gatk_files_dir+"bamout/"+specie+"/"
    bamout = bamout_dir+chrom+"."+start+"."+end+".bam"
    read_group = read_group_dir+individual+".rg.txt"

    bamout_individual_dir = bamout_files_dir+individual+"/"

    if not os.path.exists(bamout_individual_dir):
      os.makedirs(bamout_individual_dir)

    bamout_individual_bam = bamout_individual_dir+chrom+"."+start+"."+end+".bam"
    bamout_individual_bai = bamout_individual_dir+chrom+"."+start+"."+end+".bai"
    completed = compl_folder+"SplitBamOutRegion_"+individual+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

    prev_completed = compl_folder+"GetReadGroup_"+individual+".COMPLETED"

    gwf.target_from_template("SplitBamOutRegion_"+individual+"_"+chrom+"_"+start+"_"+end,
            split_bamout_region(read_group = read_group, bamout_individual_bam = bamout_individual_bam, bamout_individual_bai = bamout_individual_bai, bamout=bamout, prev_completed = prev_completed, completed = completed))

for individual in set(itertools.chain(*GIBBON_FAMILIES)):
  for chrom, start, end in GIBBON_REGIONS:

    chrom = str(chrom)
    start = str(start)
    end = str(end)

    specie = "gibbon"
    bamout_dir = new_gatk_files_dir+"bamout/"+specie+"/"
    bamout = bamout_dir+chrom+"."+start+"."+end+".bam"
    read_group = read_group_dir+individual+".rg.txt"

    bamout_individual_dir = bamout_files_dir+individual+"/"

    if not os.path.exists(bamout_individual_dir):
      os.makedirs(bamout_individual_dir)

    bamout_individual_bam = bamout_individual_dir+chrom+"."+start+"."+end+".bam"
    bamout_individual_bai = bamout_individual_dir+chrom+"."+start+"."+end+".bai"
    completed = compl_folder+"SplitBamOutRegion_"+individual+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

    prev_completed = compl_folder+"GetReadGroup_"+individual+".COMPLETED"

    gwf.target_from_template("SplitBamOutRegion_"+individual+"_"+chrom+"_"+start+"_"+end,
            split_bamout_region(read_group = read_group, bamout_individual_bam = bamout_individual_bam, bamout_individual_bai = bamout_individual_bai,bamout=bamout, prev_completed = prev_completed, completed = completed))

for child, father, mother in CHIMP_FAMILIES:
    for chrom, start, end in CHIMP_REGIONS:

      chrom = str(chrom)
      start = str(start)
      end = str(end)

      individual = child

      bamout_individual_dir = bamout_files_dir+individual+"/"
      bam = bamout_individual_dir+chrom+"."+start+"."+end+".bam"
      bai = bamout_individual_dir+chrom+"."+start+"."+end+".bam.bai"
      completed = compl_folder+"GetIndexBam_"+individual+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

      prev_completed = compl_folder+"SplitBamOutRegion_"+individual+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

      gwf.target_from_template("GetIndexBamRegions_"+individual+"_"+chrom+"_"+start+"_"+end,
              samtools_index_bam(bam = bam, bai = bai, prev_completed = prev_completed, completed = completed))

      filtered_bam_father = bamout_files_dir+father+"/"+chrom+"."+start+"."+end+".bam"
      filtered_bam_mother = bamout_files_dir+mother+"/"+chrom+"."+start+"."+end+".bam"
      filtered_bam_child = bamout_files_dir+child+"/"+chrom+"."+start+"."+end+".bam"

      family_coverage_child_dir = family_coverage_dir+child+"/"
      if not os.path.exists(family_coverage_child_dir):
        os.makedirs(family_coverage_child_dir)

      family_coverage = family_coverage_child_dir+chrom+"."+start+"."+end+".txt"
      completed = compl_folder+"GetFamilyCoverageRegion_"+father+"_"+mother+"_"+child+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

      prev_completed = compl_folder+"GetIndexBam_"+individual+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

      gwf.target_from_template("GetFamilyCoverageRegion_"+father+"_"+mother+"_"+child+"_"+chrom+"_"+start+"_"+end,
                get_family_coverage_region(filtered_bam_father = filtered_bam_father, filtered_bam_mother = filtered_bam_mother, filtered_bam_child = filtered_bam_child, family_coverage = family_coverage, prev_completed = prev_completed, completed = completed))

      #?? target('combine_coverage_region_' + child + '_ autosomes') << \
          #combine_coverage_region(child=child, regions=CHIMP_AUTOSOME_REGIONS, outname='autosomes')

for child, father, mother in BONOBO_FAMILIES:
    for chrom, start, end in BONOBO_REGIONS:

      chrom = str(chrom)
      start = str(start)
      end = str(end)

      individual = child

      bamout_individual_dir = bamout_files_dir+individual+"/"
      bam = bamout_individual_dir+chrom+"."+start+"."+end+".bam"
      bai = bamout_individual_dir+chrom+"."+start+"."+end+".bam.bai"
      completed = compl_folder+"GetIndexBam_"+individual+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

      prev_completed = compl_folder+"SplitBamOutRegion_"+individual+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

      gwf.target_from_template("GetIndexBamRegions_"+individual+"_"+chrom+"_"+start+"_"+end,
              samtools_index_bam(bam = bam, bai = bai, prev_completed = prev_completed, completed = completed))

for child, father, mother in GORILLA_FAMILIES:
    for chrom, start, end in GORILLA_REGIONS:

      chrom = str(chrom)
      start = str(start)
      end = str(end)

      individual = child

      bamout_individual_dir = bamout_files_dir+individual+"/"
      bam = bamout_individual_dir+chrom+"."+start+"."+end+".bam"
      bai = bamout_individual_dir+chrom+"."+start+"."+end+".bam.bai"
      completed = compl_folder+"GetIndexBam_"+individual+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

      prev_completed = compl_folder+"SplitBamOutRegion_"+individual+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

      gwf.target_from_template("GetIndexBamRegions_"+individual+"_"+chrom+"_"+start+"_"+end,
              samtools_index_bam(bam = bam, bai = bai, prev_completed = prev_completed, completed = completed))

for child, father, mother in ORANGUTAN_FAMILIES:
    for chrom in ORANGUTAN_CHROMOSOMES:

      chrom = str(chrom)
      start = str(start)
      end = str(end)

      individual = child

      bamout_individual_dir = bamout_files_dir+individual+"/"
      bam = bamout_individual_dir+chrom+".bam"
      bai = bamout_individual_dir+chrom+".bam.bai"
      completed = compl_folder+"GetIndexBam_"+individual+"_"+chrom+".COMPLETED"

      compl_folder+"SplitBamOutRegion_"+individual+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

      gwf.target_from_template("GetIndexBamChroms_"+individual+"_"+chrom,
              samtools_index_bam(bam = bam, bai = bai, prev_completed = prev_completed, completed = completed))

for child, father, mother in GIBBON_FAMILIES:
    for chrom in GIBBON_CHROMOSOMES:

      chrom = str(chrom)
      start = str(start)
      end = str(end)

      individual = child

      bamout_individual_dir = bamout_files_dir+individual+"/"
      bam = bamout_individual_dir+chrom+".bam"
      bai = bamout_individual_dir+chrom+".bam.bai"
      completed = compl_folder+"GetIndexBam_"+individual+"_"+chrom+".COMPLETED"

      compl_folder+"SplitBamOutRegion_"+individual+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

      gwf.target_from_template("GetIndexBamChroms_"+individual+"_"+chrom,
              samtools_index_bam(bam, bai, prev_completed, completed))

minQ=str(20)
minq=str(10)

for child, father, mother in CHIMP_FAMILIES:

  chromosomes = CHIMP_AUTOSOMES
  outname = "autosomes"

  new_family_coverage = family_coverage_dir+child+"_minQ"+minQ+"_minq"+minq+"_type_"+outname+".txt"
  completed = compl_folder+"CombineCoverageContextWr_"+child+"_minQ"+minQ+"_minq"+minq+"_type_"+outname+".COMPLETED"

  prev_completed = compl_folder+"GetFamilyCoverageRegion_"+father+"_"+mother+"_"+child+"_"+str(CHIMP_REGIONS[len(CHIMP_REGIONS)-1][0])+"_"+str(CHIMP_REGIONS[len(CHIMP_REGIONS)-1][1])+"_"+str(CHIMP_REGIONS[len(CHIMP_REGIONS)-1][2])+".COMPLETED"
  input_files = " ".join([family_coverage_dir+child+'/minQ'+minQ+'_minq'+minq+'_type' + '_' + chrom + '.txt' for chrom in chromosomes])

  gwf.target_from_template("CombineCoverageContextWr_"+child+"_minQ"+minQ+"_minq"+minq+"_type_"+outname,
          combine_coverage_context_wr(input_files=input_files, child=child, minQ=minQ, minq=minq, chromosomes=chromosomes, new_family_coverage=new_family_coverage, prev_completed = prev_completed, completed=completed))

  for chrom in CHIMP_AUTOSOMES:

    refname = CHIMP_REF

    filtered_bam_father = filtered_bam_dir+father+".bam"
    filtered_bai_father = filtered_bam_dir+father+".bam.bai"
    filtered_bam_mother = filtered_bam_dir+mother+".bam"
    filtered_bai_mother = filtered_bam_dir+mother+".bam.bai"
    filtered_bam_child = filtered_bam_dir+child+".bam"
    filtered_bai_child = filtered_bam_dir+child+".bam.bai"
    family_coverage_child_dir = family_coverage_dir+child+"/"
    family_coverage = family_coverage_dir+chrom+"."+start+"."+end+".txt"
    completed = compl_folder+"GetFamilyCoverageContextWr_"+child+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

    prev_completed = compl_folder+"CombineCoverageContextWr_"+child+"_minQ"+minQ+"_minq"+minq+"_type_"+outname+".COMPLETED"

    gwf.target_from_template("GetFamilyCoverageContextWr_"+child+"_"+chrom+"_"+start+"_"+end,
            get_family_coverage_context_wr(filtered_bam_father, filtered_bai_father, filtered_bam_mother, filtered_bai_mother, filtered_bam_child, filtered_bai_child, family_coverage, minQ, minq, refname, chrom, prev_completed, completed))

for child, father, mother in BONOBO_FAMILIES:

  chromosomes = BONOBO_AUTOSOMES
  outname = "autosomes"

  input_list = ""
  new_family_coverage = family_coverage_dir+child+"_minQ"+minQ+"_minq"+minq+"_type_"+outname+".txt"
  completed = compl_folder+"CombineCoverageContextWr_"+child+"_minQ"+minQ+"_minq"+minq+"_type_"+outname+".COMPLETED"

  prev_completed = compl_folder+"GetFamilyCoverageRegion_"+father+"_"+mother+"_"+child+"_"+str(BONOBO_REGIONS[len(BONOBO_REGIONS)-1][0])+"_"+str(BONOBO_REGIONS[len(BONOBO_REGIONS)-1][1])+"_"+str(BONOBO_REGIONS[len(BONOBO_REGIONS)-1][2])+".COMPLETED"
  input_files = " ".join([family_coverage_dir+child+'/minQ'+minQ+'_minq'+minq+'_type' + '_' + chrom + '.txt' for chrom in chromosomes])

  gwf.target_from_template("CombineCoverageContextWr_"+child+"_minQ"+minQ+"_minq"+minq+"_type_"+outname,
          combine_coverage_context_wr(input_list=input_list, child=child, minQ=minQ, minq=minq, chromosomes=chromosomes, new_family_coverage=new_family_coverage, prev_completed=prev_completed, completed=completed))

  for chrom in BONOBO_AUTOSOMES:

    refname = BONOBO_REF

    bamout_father_dir = bamout_files_dir+father+"/"
    filtered_bam_father = bamout_father_dir+chrom+"."+start+"."+end+".bam"
    filtered_bai_father = bamout_father_dir+chrom+"."+start+"."+end+".bam.bai"
    bamout_mother_dir = bamout_files_dir+mother+"/"
    filtered_bam_mother = bamout_mother_dir+chrom+"."+start+"."+end+".bam"
    filtered_bai_mother = bamout_mother_dir+chrom+"."+start+"."+end+".bam.bai"
    bamout_child_dir = bamout_files_dir+child+"/"
    filtered_bam_child = bamout_child_dir+chrom+"."+start+"."+end+".bam"
    filtered_bai_child = bamout_child_dir+chrom+"."+start+"."+end+".bam.bai"
    family_coverage_child_dir = family_coverage_dir+child+"/"
    family_coverage = family_coverage_dir+chrom+"."+start+"."+end+".txt"
    completed = compl_folder+"GetFamilyCoverageContextWr_"+child+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

    prev_completed = compl_folder+"CombineCoverageContextWr_"+child+"_minQ"+minQ+"_minq"+minq+"_type_"+outname+".COMPLETED"

    gwf.target_from_template("GetFamilyCoverageContextWr_"+child+"_"+chrom+"_"+start+"_"+end,
            get_family_coverage_context_wr(filtered_bam_father, filtered_bai_father, filtered_bam_mother, filtered_bai_mother, filtered_bam_child, filtered_bai_child, family_coverage, minQ, minq, refname, chrom, prev_completed, completed))

for child, father, mother in GORILLA_FAMILIES:
    
  chromosomes = GORILLA_AUTOSOMES
  outname = "autosomes"

  input_list = ""
  new_family_coverage = family_coverage_dir+child+"_minQ"+minQ+"_minq"+minq+"_type_"+outname+".txt"
  completed = compl_folder+"CombineCoverageContextWr_"+child+"_minQ"+minQ+"_minq"+minq+"_type_"+outname+".COMPLETED"

  prev_completed = compl_folder+"GetFamilyCoverageRegion_"+father+"_"+mother+"_"+child+"_"+str(GORILLA_REGIONS[len(GORILLA_REGIONS)-1][0])+"_"+str(GORILLA_REGIONS[len(GORILLA_REGIONS)-1][1])+"_"+str(GORILLA_REGIONS[len(GORILLA_REGIONS)-1][2])+".COMPLETED"
  input_files = " ".join([family_coverage_dir+child+'/minQ'+minQ+'_minq'+minq+'_type' + '_' + chrom + '.txt' for chrom in chromosomes])

  gwf.target_from_template("CombineCoverageContextWr_"+child+"_minQ"+minQ+"_minq"+minq+"_type_"+outname,
          combine_coverage_context_wr(input_list=input_list, child=child, minQ=minQ, minq=minq, chromosomes=chromosomes, new_family_coverage=new_family_coverage, prev_completed=prev_completed, completed=completed))

  for chrom in GORILLA_AUTOSOMES:

    refname = GORILLA_REF

    bamout_father_dir = bamout_files_dir+father+"/"
    filtered_bam_father = bamout_father_dir+chrom+"."+start+"."+end+".bam"
    filtered_bai_father = bamout_father_dir+chrom+"."+start+"."+end+".bam.bai"
    bamout_mother_dir = bamout_files_dir+mother+"/"
    filtered_bam_mother = bamout_mother_dir+chrom+"."+start+"."+end+".bam"
    filtered_bai_mother = bamout_mother_dir+chrom+"."+start+"."+end+".bam.bai"
    bamout_child_dir = bamout_files_dir+child+"/"
    filtered_bam_child = bamout_child_dir+chrom+"."+start+"."+end+".bam"
    filtered_bai_child = bamout_child_dir+chrom+"."+start+"."+end+".bam.bai"
    family_coverage_child_dir = family_coverage_dir+child+"/"
    family_coverage = family_coverage_dir+chrom+"."+start+"."+end+".txt"
    completed = compl_folder+"GetFamilyCoverageContextWr_"+child+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

    prev_completed = compl_folder+"CombineCoverageContextWr_"+child+"_minQ"+minQ+"_minq"+minq+"_type_"+outname+".COMPLETED"

    gwf.target_from_template("GetFamilyCoverageContextWr_"+child+"_"+chrom+"_"+start+"_"+end,
            get_family_coverage_context_wr(filtered_bam_father, filtered_bai_father, filtered_bam_mother, filtered_bai_mother, filtered_bam_child, filtered_bai_child, family_coverage, minQ, minq, refname, chrom, prev_completed, completed))

for child, father, mother in ORANGUTAN_FAMILIES:

  chromosomes = ORANGUTAN_AUTOSOMES
  outname = "autosomes"

  input_list = ""
  new_family_coverage = family_coverage_dir+child+"_minQ"+minQ+"_minq"+minq+"_type_"+outname+".txt"
  completed = compl_folder+"CombineCoverageContextWr_"+child+"_minQ"+minQ+"_minq"+minq+"_type_"+outname+".COMPLETED"

  prev_completed = compl_folder+"GetFamilyCoverageRegion_"+father+"_"+mother+"_"+child+"_"+str(ORANGUTAN_REGIONS[len(ORANGUTAN_REGIONS)-1][0])+"_"+str(ORANGUTAN_REGIONS[len(ORANGUTAN_REGIONS)-1][1])+"_"+str(ORANGUTAN_REGIONS[len(ORANGUTAN_REGIONS)-1][2])+".COMPLETED"
  input_files = " ".join([family_coverage_dir+child+'/minQ'+minQ+'_minq'+minq+'_type' + '_' + chrom + '.txt' for chrom in chromosomes])

  gwf.target_from_template("CombineCoverageContextWr_"+child+"_minQ"+minQ+"_minq"+minq+"_type_"+outname,
          combine_coverage_context_wr(input_list=input_list, child=child, minQ=minQ, minq=minq, chromosomes=chromosomes, new_family_coverage=new_family_coverage, prev_completed = prev_completed, completed=completed))

  for chrom in ORANGUTAN_AUTOSOMES:

    refname = ORANGUTAN_REF

    bamout_father_dir = bamout_files_dir+father+"/"
    filtered_bam_father = bamout_father_dir+chrom+"."+start+"."+end+".bam"
    filtered_bai_father = bamout_father_dir+chrom+"."+start+"."+end+".bam.bai"
    bamout_mother_dir= bamout_files_dir+mother+"/"
    filtered_bam_mother = bamout_mother_dir+chrom+"."+start+"."+end+".bam"
    filtered_bai_mother = bamout_mother_dir+chrom+"."+start+"."+end+".bam.bai"
    bamout_child_dir = bamout_files_dir+child+"/"
    filtered_bam_child = bamout_child_dir+chrom+"."+start+"."+end+".bam"
    filtered_bai_child = bamout_child_dir+chrom+"."+start+"."+end+".bam.bai"
    family_coverage_child_dir = family_coverage_dir+child+"/"
    family_coverage = family_coverage_dir+chrom+"."+start+"."+end+".txt"
    completed = compl_folder+"GetFamilyCoverageContextWr_"+child+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

    prev_completed = compl_folder+"CombineCoverageContextWr_"+child+"_minQ"+minQ+"_minq"+minq+"_type_"+outname+".COMPLETED"

    gwf.target_from_template("GetFamilyCoverageContextWr_"+child+"_"+chrom+"_"+start+"_"+end,
            get_family_coverage_context_wr(filtered_bam_father, filtered_bai_father, filtered_bam_mother, filtered_bai_mother, filtered_bam_child, filtered_bai_child, family_coverage, minQ, minq, refname, chrom, prev_completed, completed))

for child, father, mother in GIBBON_FAMILIES:

  chromosomes = GIBBON_AUTOSOMES
  outname = "autosomes"

  input_list = ""
  new_family_coverage = family_coverage_dir+child+"_minQ"+minQ+"_minq"+minq+"_type_"+outname+".txt"
  completed = compl_folder+"CombineCoverageContextWr_"+child+"_minQ"+minQ+"_minq"+minq+"_type_"+outname+".COMPLETED"

  prev_completed = compl_folder+"GetFamilyCoverageRegion_"+father+"_"+mother+"_"+child+"_"+str(GIBBON_REGIONS[len(GIBBON_REGIONS)-1][0])+"_"+str(GIBBON_REGIONS[len(GIBBON_REGIONS)-1][1])+"_"+str(GIBBON_REGIONS[len(GIBBON_REGIONS)-1][2])+".COMPLETED"
  input_files = " ".join([family_coverage_dir+child+'/minQ'+minQ+'_minq'+minq+'_type' + '_' + chrom + '.txt' for chrom in chromosomes])

  gwf.target_from_template("CombineCoverageContextWr_"+child+"_minQ"+minQ+"_minq"+minq+"_type_"+outname,
          combine_coverage_context_wr(input_list=input_list, child=child, minQ=minQ, minq=minq, chromosomes=chromosomes, new_family_coverage=new_family_coverage, prev_completed=prev_completed, completed=completed))

  for chrom in GIBBON_AUTOSOMES:

    refname = GIBBON_REF

    bamout_father_dir = bamout_files_dir+father+"/"
    filtered_bam_father = bamout_father_dir+chrom+"."+start+"."+end+".bam"
    filtered_bai_father = bamout_father_dir+chrom+"."+start+"."+end+".bam.bai"
    bamout_mother_dir = bamout_files_dir+mother+"/"
    filtered_bam_mother = bamout_mother_dir+chrom+"."+start+"."+end+".bam"
    filtered_bai_mother = bamout_mother_dir+chrom+"."+start+"."+end+".bam.bai"
    bamout_child_dir = bamout_files_dir+child+"/"
    filtered_bam_child = bamout_child_dir+chrom+"."+start+"."+end+".bam"
    filtered_bai_child = bamout_child_dir+chrom+"."+start+"."+end+".bam.bai"
    family_coverage_child_dir = family_coverage_dir+child+"/"
    family_coverage = family_coverage_dir+chrom+"."+start+"."+end+".txt"
    completed = compl_folder+"GetFamilyCoverageContextWr_"+child+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

    prev_completed = compl_folder+"CombineCoverageContextWr_"+child+"_minQ"+minQ+"_minq"+minq+"_type_"+outname+".COMPLETED"

    gwf.target_from_template("GetFamilyCoverageContextWr_"+child+"_"+chrom+"_"+start+"_"+end,
            get_family_coverage_context_wr(filtered_bam_father, filtered_bai_father, filtered_bam_mother, filtered_bai_mother, filtered_bam_child, filtered_bai_child, family_coverage, minQ, minq, refname, chrom, prev_completed, completed))

