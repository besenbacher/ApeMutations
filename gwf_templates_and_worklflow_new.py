######################################################################################
#A.SET ENVIRONMENT AND IMPORT SYS,OS
###################################################################################### 
from gwf import Workflow
import sys, os
gwf = Workflow()

######################################################################################
#B.FUNCTIONS AND TEMPLATES
######################################################################################

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
  output = [ref_genome_dict,ref_genome_fai,completed]
  options = {"memory": "1g","walltime":"01:00:00"}
  spec='''
  source activate prjMutRates
  source /com/extra/samtools/LATEST/load.sh
  
  samtools faidx ref-genomes/{ref_genome_fasta}

  source /com/extra/java/8/load.sh
  source /com/extra/picard/LATEST/load.sh
  
  picard CreateSequenceDictionary R=ref-genomes/{ref_genome_fasta} O=ref-genomes/{ref_genome_dict}

  echo "Completed at "$(date) > {completed}'''.format(ref_genome_fasta=ref_genome_fasta,ref_genome_dict=ref_genome_dict, completed=completed)
  
  return inputs, outputs, options, spec

def samtools_index_bam(bam,bai,completed):
  inputs = [bam]
  output = [bai,completed]
  options = {"memory": "1g","walltime":"00:59:00"}
  spec='''
  source activate prjMutRates
  source /com/extra/samtools/LATEST/load.sh
  
  samtools index {bam}

  echo "Completed at "$(date) > {completed}'''.format(bam=bam, completed=completed)
  
  return inputs, outputs, options, spec

def collect_realign_regions(ref_genome_fasta, ref_genome_dict, merged_bam, merged_bai, realign_intervals, completed):
  inputs = [ref_genome_fasta, ref_genome_dict, merged_bam, merged_bai]
  output = [realign_intervals, completed]
  options = {"cores": 16}
  spec='''
  source activate prjMutRates
  source /com/extra/java/8/load.sh
  source /com/extra/GATK/3.8/load.sh

  gatk -nt 16 -T RealignerTargetCreator -R ref-genomes/{ref_genome_fasta} \
     -I merged_bam_files/{merged_bam_files} \
     -o new_gatk_files/{realign_intervals}

  echo "Completed at "$(date) > {completed}'''.format(ref_genome_fasta=ref_genome_fasta, merged_bam=merged_bam, realign_intervals=realign_intervals, completed=completed)
  
  return inputs, outputs, options, spec

def realign(ref_genome_fasta, ref_genome_dict, merged_bam, merged_bai, realign_intervals, realign_bam, completed):
  inputs = [ref_genome_fasta, ref_genome_dict, merged_bam, merged_bai, realign_intervals]
  output = [realign_bam, realign_bai, completed]
  options = {"cores": 16}
  spec='''
  source /com/extra/java/8/load.sh
  source /com/extra/GATK/3.8/load.sh

  gatk -T IndelRealigner -R ref-genomes/{ref_genome_fasta} \
     --filter_bases_not_stored \
     -I merged_bam_files/{merged_bam} \
     -targetIntervals new_gatk_files/{realign_intervals} \
     -o new_gatk_files/{realign_bam}

  echo "Completed at "$(date) > {completed}'''.format(ref_genome_fasta=ref_genome_fasta, merged_bam=merged_bam, realign_intervals=realign_intervals, realign_bam = realign_bam, completed=completed)
  
  return inputs, outputs, options, spec

def calc_recalibrate_info(ref_genome_fasta, ref_genome_dict, realign_bam, realign_bai, knownlist, recalibration_report, completed):
  inputs = [ref_genome_fasta, ref_genome_dict, realign_bam, realign_bai, knownlist]
  output = [recalibration_report, completed]
  options = {"cores": 16}
  spec='''
  source /com/extra/java/8/load.sh
  source /com/extra/GATK/3.8/load.sh

  gatk -T BaseRecalibrator \
     -nct 16 \
     -R {ref_genome_fasta} \
     {knownlist}  \
     -I {realign_bam} \
     -o {recalibration_report}
  
  echo "Completed at "$(date) > {completed}'''.format(ref_genome_fasta=ref_genome_fasta, realign_bam = realign_bam, knownlist=knownlist, recalibration_report=recalibration_report, completed=completed)
  
  return inputs, outputs, options, spec

def call_calc_recalibrate_info(**arguments):
  arguments['knownlist'] = ' '.join('-knownSites ' + x for x in arguments['knownlist'])
  return calc_recalibrate_info(**arguments)

def recalibrate(ref_genome_fasta, ref_genome_dict, realign_bam, realign_bai, recalibration_report, recalibrated_bam, recalibrated_bai, completed):
  inputs = [ref_genome_fasta, ref_genome_dict, realign_bam, realign_bai, recalibration_report]
  outputs = [recalibrated_bam, recalibrated_bai, completed]
  options = {"cores": 16, "memory":"50g"}
  spec='''
  source /com/extra/java/8/load.sh
  source /com/extra/GATK/3.8/load.sh

  java -Djava.io.tmpdir=/scratch/$PBS_JOBID/tmp \
     -Djava.awt.headless=true \
     -Xmx50g \
     -jar /com/extra/GATK/3.8/jar-bin/GenomeAnalysisTK.jar \
     -T PrintReads \
     -nct 16 \
     -R {ref_genome_fasta} \
     -I {realign_bam} \
     -BQSR {recalibration_report} \
     -o {recalibrated_bam}

  echo "Completed at "$(date) > {completed}'''.format(ref_genome_fasta=ref_genome_fasta, ref_genome_dict=ref_genome_dict, realign_bam=realign_bam, realign_bai=realign_bai, recalibration_report=recalibration_report, recalibrated_bam=recalibrated_bam, recalibrated_bai=recalibrated_bai, completed=completed)

  return inputs, outputs, options, spec

def filter_bam_files(recalibrated_bam, recalibrated_bai, ref_genome_fasta, filtered_bam, completed):
  inputs = [recalibrated_bam, recalibrated_bai, ref_genome_fasta]
  outputs = [filtered_bam, completed]
  options = {"cores": 16, "walltime":"12:00:00"}
  spec='''
  source /com/extra/java/8/load.sh
  source /com/extra/GATK/3.8/load.sh

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

def haplotype_caller_region_bo(ref_genome_fasta, haplocaller_vcf, haplocaller_vcf_idx, bamout, individuals, chr, start, end, includelist, completed, extra_options):
  inputs = [ref_genome_fasta]
  output = [haplocaller_vcf, haplocaller_vcf_idx, bamout, completed]
  options = {"cores": 1, "memory": "16g","walltime":"160:00:00"}
  spec='''
  source /com/extra/java/8/load.sh
  source /com/extra/GATK/3.8/load.sh  

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

  echo "Completed at "$(date) > {completed}'''.format(ref_genome_fasta=ref_genome_fasta, haplocaller_vcf=haplocaller_vcf, haplocaller_vcf_idx=haplocaller_vcf_idx, bamout=bamout, completed=completed)

  return inputs, outputs, options, spec

def call_haplotype_caller_region_bo(**arguments):
  assert 'individuals' in arguments
  bamfiles = [new_gatk_files_dir+individual+'.recalibrated.bam' for individual in arguments['individuals']]
  extra_options = {'input':bamfiles}
  arguments['includelist'] = ' '.join('-I ' + x for x in bamfiles)
  (options, spec) = haplotype_caller_region_bo(**arguments)
  return (add_options(options,extra_options), spec)

def combine_regions(split_file, species_haplocaller_vcf, completed):
  inputs = [split_file]
  output = [species_haplocaller_vcf, completed]
  options = {"walltime":"10:00:00"}
  spec='''
  source /com/extra/vcftools/LATEST/load.sh

  vcf-concat $(cat {split_file} | fgrep -v chrX | fgrep -v chrUn | fgrep -v chrM | fgrep -v chrY | fgrep -v hap | ~/Scripts/gorsort.sh | awk -vORS=" " '{{print "new_vcf_files_dir/{specie}/haplocaller.raw."$1"."$2"."$3".vcf"}}') > {species_haplocaller_vcf}

  echo "Completed at "$(date) > {completed}'''.format(split_file=split_file, species_haplocaller_vcf=species_haplocaller_vcf, completed=completed)

  return inputs, outputs, options, spec

def call_combine_regions(**arguments):
    vcf_files = [new_vcf_files_dir + arguments["specie"] + '/haplocaller.raw.' + chrom + '.' + str(start) + '.' + str(end) + '.vcf' for (chrom,start,end) in arguments['regions']]
    extra_options = {'input':vcf_files}
    (options, spec) = combine_regions(**arguments)
    return (add_options(options, extra_options), spec)

refname=CHIMP_REF, specie="chimpanzees", regions=CHIMP_REGIONS

def get_readgroup(merged_bam_rg, read_group, completed):
  inputs = [merged_bam_rg]
  ouptut = [read_group, completed]
  options = {"memory":"1g","walltime":"00:59:00"}
  spec='''

  cat {merged_bam_rg} | cut -f2 | cut -d: -f2 > {read_group}

  echo "Completed at "$(date) > {completed}'''.format(merged_bam_rg=merged_bam_rg, read_groups=read_groups, completed=completed)

  return inputs, outputs, options, spec

def split_bamout_region(read_group, bamout, bamout_individual_bam, bamout_individual_bai, completed):
  inputs = [read_groups, bamout]
  ouptut = [bamout_individual_bam, bamout_individual_bai, completed]
  options = {"memory":"4g","walltime":"10:00:00"}
  spec='''
  source /com/extra/samtools/LATEST/load.sh

  samtools view -R {read_group} -b -o {bamout_individual_bam} {bamout}
  samtools index {bamout_individual_bai}

  echo "Completed at "$(date) > {completed}'''.format(read_group = read_group, bamout = bamout, bamout_individual_bam = bamout_individual_bam, bamout_individual_bai = bamout_individual_bai, completed)

  return inputs, outputs, options, spec

#inkluderer info om hvorvidt basen er repeat_masket

def get_family_coverage_region(filtered_bam_father, filtered_bam_mother, filtered_bam_child, family_coverage, completed):
  inputs = [filtered_bam_father, filtered_bam_mother, filtered_bam_child]
  output = [family_coverage, completed]
  options = {"memory":"1g","walltime":"10:00:00"}
  spec='''
  source /com/extra/samtools/LATEST/load.sh

  samtools depth {filtered_bam_father} {filtered_bam_mother} {filtered_bam_child} | awk '$3>=5 && $3<=150 && $4>=5 && $4<=150 && $5>=5 && $5<=150' | cut -f3- | ~/Scripts/table.py > {family_coverage}

  echo "Completed at "$(date) > {completed}'''.format(filtered_bam_father = filtered_bam_father, filtered_bam_mother = filtered_bam_mother, filter_bam_child = filtered_bam_child, completed = completed)

  return inputs, outputs, options, spec

def combine_coverage_context_wr(input_list, new_family_coverage, completed):
  inputs = [input_list]
  output = [new_family_coverage, completed]
  options = {"memory":"8g","walltime":"10:00:00"}
  spec='''

  cat {input_list} | ~/Scripts/combine_tables.py > {new_family_coverage}

  echo "Completed at "$(date) > {completed}'''.format(input_list = input_list, new_family_coverage = new_family_coverage, completed = completed)

  return inputs, outputs, options, spec

def call_combine_coverage_context_wr(**arguments):
    chromosomes = arguments['chromosomes']
    input_files = ['new_family_coverage_wr/{child}/minQ{minQ}_minq{minq}_type'.format(**arguments) + '_' + chrom + '.txt' for chrom in chromosomes]
    arguments['input_list'] = ' '.join(input_files)
    extra_options = {'input':input_files}
    (options, spec) = _combine_coverage_context_wr(**arguments)
    return (add_options(options, extra_options), spec)

get_family_coverage_chrom = \
    template(input=['bamout_files/{father}/{chrom}.bam',
                    'bamout_files/{mother}/{chrom}.bam',
                    'bamout_files/{child}/{chrom}.bam'],
             output=['new_family_coverage/{child}/{chrom}.txt'],
             account="DanishPanGenome") << '''
source /com/extra/samtools/LATEST/load.sh

mkdir -p new_family_coverage/{child}/

samtools depth bamout_files/{father}/{chrom}.bam bamout_files/{mother}/{chrom}.bam bamout_files/{child}/{chrom}.bam | awk '$3>=5 && $3<=150 && $4>=5 && $4<=150 && $5>=5 && $5<=150' | cut -f3- | ~/Scripts/table.py > new_family_coverage/{child}/{chrom}.txt
'''

def get_family_coverage_context_wr(filtered_bam_father, filtered_bai_father, filtered_bam_mother, filtered_bai_mother, filtered_bam_child, filtered_bai_child, family_coverage, minQ, minq, refname, completed):
  inputs = [filtered_bam_father, filtered_bai_father, filtered_bam_mother, filtered_bai_mother, filtered_bam_child, filtered_bai_child]
  output = [family_coverage, completed]
  options = {"walltime":"48:00:00"}
  spec='''
  source /com/extra/samtools/LATEST/load.sh
  source /com/extra/python/LATEST/load.sh

  samtools depth -Q{minQ} -q{minq} {filtered_bam_father} {filtered_bam_mother} {filtered_bam_child} -r{chrom}| awk '$3>=5 && $3<=200 && $4>=5 && $4<=200 && $5>=5 && $5<=200' | ./python_scripts/add_context_type_w_repeat.py {refname} | cut -f3- | ~/Scripts/table.py > {family_coverage}

  echo "Completed at "$(date) > {completed}'''.format(filtered_bam_father = filtered_bam_father, filter_bai_father = filtered_bam_father, filtered_bam_mother = filtered_bam_mother, filter_bai_mother = filtered_bam_mother, filtered_bam_child = filtered_bam_child, filter_bai_child = filtered_bam_child, family_coverage = family_coverage, minQ=minQ, minq=minq, refname = refname, completed = completed)

  return inputs, outputs, options, spec

######################################################################################
#C.PATHS AND VARIABLES
###################################################################################### 

pwd                 = 
out                 = pwd+"out/"
compl_folder        = pwd+"out/COMPLETED/"

ref_genomes_dir     = pwd+"ref_genomes/"
merged_bams_dir     = pwd+"merged_bams/"
filtered_bam_dir    = pwd+"filtered_bams/" 
new_gatk_files_dir  = pwd+"new_gatk_files/" 
split_files_dir     = pwd+"split_files/"
new_vcf_files_dir   = pwd+"new_vcf_files/"
read_group_dir      = pwd+"read_groups/"
bamout_files_dir    = pwd+"bamout_files/"
family_coverage_dir = pwd+"new_family_coverage_wr"

ORANGUTAN_REF = 'ponAbe2'
CHIMP_REF     = 'panTro5' 
GORILLA_REF   = 'gorGor4'

ORANGUTANS = ['Buschi', 'Moni', 'Masala', 'Farida', 'Ogan',  'Aris', 'Schubbi', 'Pongo']
CHIMPS     = ['Frits', 'Carolina', 'Carl', 'Simliki', 'ERR466113','ERR466114','ERR466115','ERR466116','ERR466117','ERR466118','ERR466119','ERR466120','ERR466121']
GORILLAS   = ['Banjo', 'Mimi', 'Mutasi', 'Mawenzi', 'Efata', 'Undi']

CHIMP_FAMILIES = read_family_description('family_description/chimpanzees.txt')
GORILLA_FAMILIES = read_family_description('family_description/gorillas.txt')
ORANGUTAN_FAMILIES = read_family_description('family_description/orangutans.txt')

#Order: (Child, Father, Mother)

ALL_APES = ORANGUTANS + CHIMPS + GORILLAS


orang_known = ['../KnownVariants/abelii.2012-08-17.snps.dropsamples.autos.vqsr99.vcf',
               '../KnownVariants/abelii.2012-08-17.snps.dropsamples.chrX.vcf',
               '../KnownVariants/pygmaeus.2012-08-17.snps.dropsamples.autos.vqsr99.vcf',
               '../KnownVariants/pygmaeus.2012-08-17.snps.dropsamples.chrX.vcf']

chimp_known = ['../faststorage/KnownVariants/pantro.snps.2012-11-25.autos.99vqsr.beagle_' + CHIMP_REF +'_lift.vcf',
               '../faststorage/KnownVariants/gatk.allchimpanzee.2012-05-30.snps.chrX.vqsr99_' + CHIMP_REF + '_lift.vcf']
gorilla_known = ['../faststorage/KnownVariants/emitall-var.dropsamples.autos.vqsr99.2012-04-26.abfilter_' + GORILLA_REF  + '_lift_sorted.vcf']


######################################################################################
#D.CODE
###################################################################################### 

#First we need to index the reference genomes with Picard
#refGenomes should be in ./ref-genomes

for ref in ORANGUTAN_REF, CHIMP_REF, GORILLA_REF:

  ref_genome_fasta = ref_genomes_dir+ref+".fa"
  ref_genome_fai = ref_genomes_dir+ref+".fa.fai"
  ref_genome_dict = ref_genomes_dir+ref+".dict"
  completed = compl_folder+"GetIndexRef_"+ref+".COMPLETED"

  gwf.target_from_template("GetIndexRef_"+ref,
            picard_index_reference(ref_genome_fasta,ref_genome_dict,ref_genome_fai,completed))

# Use GATK to process bam files before calling.

#To work with the bam files in GATK we first need them indexed...
# bamfiles should be in ./merged-bams

for individual in ALL_APES:

  merged_bam = merged_bams_dir+individual+".bam"
  merged_bai = merged_bams_dir+individual+".bam.bai"
  completed = compl_folder+"GetIndexMergedBam_"+individual+".COMPLETED"

  gwf.target_from_template("GetIndexMergedBam_"+individual,
        samtools_index_bam(merged_bam,merged_bai,completed))

  filtered_bam = filtered_bam_dir+individual+".bam"
  filtered_bai = filtered_bam_dir+individual+".bam.bai"
  completed = compl_folder+"GetIndexFilteredBam_"+individual+".COMPLETED"

  gwf.target_from_template("GetIndexFilteredBam_"+individual,
        samtools_index_bam(filtered_bam,filtered_bai,completed))


#Next, identify regions that needs to be locally re-aligned.
# intermediary output files from gatk will be put in ./gatk_files/

#Realign the bam files:

for individual in CHIMPS:

  ref_genome_fasta = ref_genomes_dir+CHIMP_REF+".fa"
  ref_genome_dict = ref_genomes_dir+CHIMP_REF+".dict"
  merged_bam = merged_bams_dir+individual+".bam"
  merged_bai = merged_bams_dir+individual+".bai"
  realign_intervals = new_gatk_files_dir+individual+".realign.intervals"
  completed = compl_folder+"CollectRealignRegions_"+individual+".COMPLETED"

  gwf.target_from_template("CollectRealignRegions_"+individual,
        collect_realign_regions(ref_genome_fasta, ref_genome_dict, merged_bam, merged_bai, realign_intervals, completed))

  realign_bam = new_gatk_files_dir+individual+".realigned.bam"
  realign_bai = new_gatk_files_dir+individual+".realigned.bai"
  completed = compl_folder+"Realign_"+individual+".COMPLETED"

  gwf.target_from_template("Realign_"+individual,
        realign(ref_genome_fasta, ref_genome_dict, merged_bam, merged_bai, realign_intervals, realign_bam, completed))

  recalibration_report = new_gatk_files_dir+individual+".recalibration_report.grp"
  completed = compl_folder+"CalcRecInfo_"+individual+".COMPLETED"

  knownlist = ""

  gwf.target_from_template("CalcRecInfo_"+individual,
        call_calc_recalibrate_info(ref_genome_fasta=ref_genome_fasta, ref_genome_dict=ref_genome_dict, realign_bam=realign_bam, realign_bai=realign_bai, knownlist=knownlist, chimp_known=chimp_known, recalibration_report=recalibration_report, completed=completed))

  recalibrated_bam = new_gatk_files_dir+individual+".recalibrated.bam"
  completed = compl_folder+"GetRecalBam_"+individual+".COMPLETED"

  gwf.target_from_template("GetRecalBam_"+individual,
        recalibrate(ref_genome_fasta, ref_genome_dict, realign_bam, realign_bai, recalibration_report, recalibrated_bam, recalibrated_bai, completed))

  filtered_bam = filtered_bam_dir+individual+".bam"
  completed = compl_folder+"FilterBam_"+individual+".COMPLETED"

  gwf.target_from_template("FilterBam_"+individual,
        filter_bam_files(recalibrated_bam, recalibrated_bai, ref_genome_fasta, filtered_bam, completed))

for individual in GORILLAS:

  ref_genome_fasta = ref_genomes_dir+GORILLA_REF+".fa"
  ref_genome_dict = ref_genomes_dir+GORILLA_REF+".dict"
  merged_bam = merged_bams_dir+individual+".bam"
  merged_bai = merged_bams_dir+individual+".bai"
  realign_intervals = new_gatk_files_dir+individual+".realign.intervals"
  completed = compl_folder+"CollectRealignRegions_"+individual+".COMPLETED"

  gwf.target_from_template("CollectRealignRegions_"+individual,
        collect_realign_regions(ref_genome_fasta, ref_genome_dict, merged_bam, merged_bai, realign_intervals, completed))

  realign_bam = new_gatk_files_dir+individual+".realigned.bam"
  realign_bai = new_gatk_files_dir+individual+".realigned.bai"
  completed = compl_folder+"Realign_"+individual+".COMPLETED"

  gwf.target_from_template("Realign_"+individual,
        realign(ref_genome_fasta, ref_genome_dict, merged_bam, merged_bai, realign_intervals, realign_bam, completed))

  recalibration_report = new_gatk_files_dir+individual+".recalibration_report.grp"
  completed = compl_folder+"CalcRecInfo_"+individual+".COMPLETED"

  knownlist = ""

  gwf.target_from_template("CalcRecInfo_"+individual,
        call_calc_recalibrate_info(ref_genome_fasta=ref_genome_fasta, ref_genome_dict=ref_genome_dict, realign_bam=realign_bam, realign_bai=realign_bai, knownlist=knownlist, chimp_known=chimp_known, recalibration_report=recalibration_report, completed=completed))

  recalibrated_bam = new_gatk_files_dir+individual+".recalibrated.bam"
  completed = compl_folder+"GetRecalBam_"+individual+".COMPLETED"

  gwf.target_from_template("GetRecalBam_"+individual,
        recalibrate(ref_genome_fasta, ref_genome_dict, realign_bam, realign_bai, recalibration_report, recalibrated_bam, recalibrated_bai, completed))

  filtered_bam = filtered_bam_dir+individual+".bam"
  completed = compl_folder+"FilterBam_"+individual+".COMPLETED"

  gwf.target_from_template("FilterBam_"+individual,
        filter_bam_files(recalibrated_bam, recalibrated_bai, ref_genome_fasta, filtered_bam, completed))

for individual in ORANGUTANS:

  ref_genome_fasta = ref_genomes_dir+ORANGUTAN_REF+".fa"
  ref_genome_dict = ref_genomes_dir+ORANGUTAN_REF+".dict"
  merged_bam = merged_bams_dir+individual+".bam"
  merged_bai = merged_bams_dir+individual+".bai"
  realign_intervals = new_gatk_files_dir+individual+".realign.intervals"
  completed = compl_folder+"CollectRealignRegions_"+individual+".COMPLETED"

  gwf.target_from_template("CollectRealignRegions_"+individual,
        collect_realign_regions(ref_genome_fasta, ref_genome_dict, merged_bam, merged_bai, realign_intervals, completed))

  realign_bam = new_gatk_files_dir+individual+".realigned.bam"
  realign_bai = new_gatk_files_dir+individual+".realigned.bai"
  completed = compl_folder+"Realign_"+individual+".COMPLETED"

  gwf.target_from_template("Realign_"+individual,
        realign(ref_genome_fasta, ref_genome_dict, merged_bam, merged_bai, realign_intervals, realign_bam, completed))

  recalibration_report = new_gatk_files_dir+individual+".recalibration_report.grp"
  completed = compl_folder+"CalcRecInfo_"+individual+".COMPLETED"

  knownlist = ""

  gwf.target_from_template("CalcRecInfo_"+individual,
        call_calc_recalibrate_info(ref_genome_fasta=ref_genome_fasta, ref_genome_dict=ref_genome_dict, realign_bam=realign_bam, realign_bai=realign_bai, knownlist=knownlist, chimp_known=chimp_known, recalibration_report=recalibration_report, completed=completed))

  recalibrated_bam = new_gatk_files_dir+individual+".recalibrated.bam"
  completed = compl_folder+"GetRecalBam_"+individual+".COMPLETED"

  gwf.target_from_template("GetRecalBam_"+individual,
        recalibrate(ref_genome_fasta, ref_genome_dict, realign_bam, realign_bai, recalibration_report, recalibrated_bam, recalibrated_bai, completed))

  filtered_bam = filtered_bam_dir+individual+".bam"
  completed = compl_folder+"FilterBam_"+individual+".COMPLETED"

  gwf.target_from_template("FilterBam_"+individual,
        filter_bam_files(recalibrated_bam, recalibrated_bai, ref_genome_fasta, filtered_bam, completed))

#Call variants using haplotype caller

def get_chromosomes(refname):
    f = open('ref_genomes_dir/'+ refname + '.fa.fai')
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

def get_regions(refname):
    f = open(split_files_dir + refname + '_split_1000_1000000.txt')
    L = []
    for line in f:
        chrom, start, end = line.split()
        L.append((chrom, int(start), int(end)))
    f.close()
    return L


CHIMP_CHROMOSOMES = get_chromosomes(CHIMP_REF)
ORANGUTAN_CHROMOSOMES = get_chromosomes(ORANGUTAN_REF)
GORILLA_CHROMOSOMES = get_chromosomes(GORILLA_REF)

ORANGUTAN_REGIONS = get_regions(ORANGUTAN_REF)
CHIMP_REGIONS = get_regions(CHIMP_REF)
GORILLA_REGIONS = get_regions(GORILLA_REF)

CHIMP_AUTOSOME_REGIONS = [(chrom, start, end) for (chrom, start, end) in CHIMP_REGIONS if not (chrom.startswith('chrX') or chrom.startswith('chrY') or chrom.startswith('chrM'))]
CHIMP_X_REGIONS = [(chrom, start, end) for (chrom, start, end) in CHIMP_REGIONS if chrom.startswith('chrX')]

GORILLA_AUTOSOME_REGIONS = [(chrom, start, end) for (chrom, start, end) in GORILLA_REGIONS if not (chrom.startswith('chrX') or chrom.startswith('chrY') or chrom.startswith('chrM'))]
GORILLA_X_REGIONS = [(chrom, start, end) for (chrom, start, end) in GORILLA_REGIONS if chrom.startswith('chrX')]

ORANGUTAN_AUTOSOME_REGIONS = [(chrom, start, end) for (chrom, start, end) in ORANGUTAN_REGIONS if not (chrom.startswith('chrX') or chrom.startswith('chrY') or chrom.startswith('chrM') or 'hap' in chrom)]
ORANGUTAN_X_REGIONS = [(chrom, start, end) for (chrom, start, end) in ORANGUTAN_REGIONS if chrom.startswith('chrX')]

ORANGUTAN_AUTOSOMES = [x for x in ORANGUTAN_CHROMOSOMES if is_autosome(x)]
GORILLA_AUTOSOMES = [x for x in GORILLA_CHROMOSOMES if is_autosome(x)]
CHIMP_AUTOSOMES = [x for x in CHIMP_CHROMOSOMES if is_autosome(x)]

for chrom, start, end in ORANGUTAN_AUTOSOME_REGIONS:

  individuals = ORANGUTANS
  ref_genome_fasta = ref_genomes_dir+ORANGUTAN_REF+".fa"
  specie = "orangutans"

  haplocaller_out_dir = new_vcf_files_dir+specie+"/"
  bamout_dir = new_gatk_files+"bamout/"+specie+"/"

  if not os.path.exists(haplocaller_out_dir):
    os.makedirs(haplocaller_out_dir)

  if not os.path.exists(bamout_dir):
    os.makedirs(bamout_dir)

  haplocaller_vcf = haplocaller_out_dir+"haplocaller.raw."+chrom+"."+start+"."+end+".vcf"
  haplocaller_vcf_idx = haplocaller_out_dir+"haplocaller.raw."+chrom+"."+start+"."+end+".vcf.idx"
  bamout = bamout_dir+chrom+"."+start+"."+end+".bam"
  completed = compl_folder+"HaplotypeCaller_"+chrom"_"+start+"_"+end+".COMPLETED"

  includelist = ""

  gwf.target_from_template("HaplotypeCaller_"+chrom"_"+start+"_"+end,
          call_haplotype_caller_region_bo(ref_genome_fasta=ref_genome_fasta, haplocaller_vcf=haplocaller_vcf, haplocaller_vcf_idx=haplocaller_vcf_idx, bamout=bamout, individuals=individuals, chrom=chrom, start=start, end=end, includelist=includelist, completed=completed))

  split_file = split_files_dir + refname + "_split_1000_1000000.txt"
  species_haplocaller_vcf = new_vcf_files_dir+specie+".haplocaller.raw.auto.vcf"
  completed = compl_folder+"CombineRegionsChimps_Gorillas.COMPLETED"

  gwf.target_from_template("CombineRegionsChimps_Gorillas",
          call_combine_regions(refname=CHIMP_REF, specie="chimpanzees", regions=CHIMP_REGIONS, split_file = split_file, species_haplocaller_vcf = species_haplocaller_vcf, completed = completed))

for chrom, start, end in GORILLA_REGIONS:

  individuals = GORILLAS
  ref_genome_fasta = ref_genomes_dir+GORILLA_REF+".fa"
  specie = "gorillas"

  haplocaller_out_dir = new_vcf_files_dir+specie+"/"
  bamout_dir = new_gatk_files+"bamout/"+specie+"/"

  if not os.path.exists(haplocaller_out_dir):
    os.makedirs(haplocaller_out_dir)

  if not os.path.exists(bamout_dir):
    os.makedirs(bamout_dir)

  haplocaller_vcf = haplocaller_out_dir+"haplocaller.raw."+chrom+"."+start+"."+end+".vcf"
  haplocaller_vcf_idx = haplocaller_out_dir+"haplocaller.raw."+chrom+"."+start+"."+end+".vcf.idx"
  bamout = bamout_dir+chrom+"."+start+"."+end+".bam"
  completed = compl_folder+"HaplotypeCaller_"+chrom"_"+start+"_"+end+".COMPLETED"

  includelist = ""

  gwf.target_from_template("HaplotypeCaller_"+chrom"_"+start+"_"+end,
          call_haplotype_caller_region_bo(ref_genome_fasta=ref_genome_fasta, haplocaller_vcf=haplocaller_vcf, haplocaller_vcf_idx=haplocaller_vcf_idx, bamout=bamout, individuals=individuals, chrom=chrom, start=start, end=end, includelist=includelist, completed=completed))

  split_file = split_files_dir + refname + "_split_1000_1000000.txt"
  species_haplocaller_vcf = new_vcf_files_dir+specie+".haplocaller.raw.auto.vcf"
  completed = compl_folder+"CombineRegionsChimps_Gorillas.COMPLETED"

  gwf.target_from_template("CombineRegionsChimps_Gorillas",
          call_combine_regions(refname=CHIMP_REF, specie="chimpanzees", regions=CHIMP_REGIONS, split_file = split_file, species_haplocaller_vcf = species_haplocaller_vcf, completed = completed))


for chrom, start, end in CHIMP_REGIONS:

  individuals = CHIMPS
  ref_genome_fasta = ref_genomes_dir+CHIMP_REF+".fa"
  specie = "chimpanzees"

  haplocaller_out_dir = new_vcf_files_dir+specie+"/"
  bamout_dir = new_gatk_files+"bamout/"+specie+"/"

  if not os.path.exists(haplocaller_out_dir):
    os.makedirs(haplocaller_out_dir)

  if not os.path.exists(bamout_dir):
    os.makedirs(bamout_dir)

  haplocaller_vcf = haplocaller_out_dir+"haplocaller.raw."+chrom+"."+start+"."+end+".vcf"
  haplocaller_vcf_idx = haplocaller_out_dir+"haplocaller.raw."+chrom+"."+start+"."+end+".vcf.idx"
  bamout = bamout_dir+chrom+"."+start+"."+end+".bam"
  completed = compl_folder+"HaplotypeCaller_"+chrom"_"+start+"_"+end+".COMPLETED"

  includelist = ""

  gwf.target_from_template("HaplotypeCaller_"+chrom"_"+start+"_"+end,
          call_haplotype_caller_region_bo(ref_genome_fasta=ref_genome_fasta, haplocaller_vcf=haplocaller_vcf, haplocaller_vcf_idx=haplocaller_vcf_idx, bamout=bamout, individuals=individuals, chrom=chrom, start=start, end=end, includelist=includelist, completed=completed))

  split_file = split_files_dir + refname + "_split_1000_1000000.txt"
  species_haplocaller_vcf = new_vcf_files_dir+specie+".haplocaller.raw.auto.vcf"
  completed = compl_folder+"CombineRegionsChimps_Chimps.COMPLETED"

  gwf.target_from_template("CombineRegionsChimps_Chimps",
          call_combine_regions(refname=CHIMP_REF, specie="chimpanzees", regions=CHIMP_REGIONS, split_file = split_file, species_haplocaller_vcf = species_haplocaller_vcf, completed = completed))

for individual in ALL_APES:

  merged_bam_rg = merged_bams_dir+individual+".rg.txt"
  read_group = read_group_dir+individual".rg.txt"
  completed = compl_folder+"GetReadGroup_"+individual+".COMPLETED"

  gwf.target_from_template("GetReadGroup_"+individual
          get_readgroup(merged_bam_rg, read_group, completed))

for individual in set(itertools.chain(*CHIMP_FAMILIES)):
  for chrom, start, end in CHIMP_REGIONS:

    specie = "chimpanzees"
    bamout_dir = new_gatk_files+"bamout/"+specie+"/"
    bamout = bamout_dir+chrom+"."+start+"."+end+".bam"
    read_group = read_group_dir+individual+".rg.txt"

    bamout_individual_dir = bamout_files_dir+individual+"/"

    if not os.path.exists(bamout_individual_dir):
      os.makedirs(bamout_individual_dir)

    bamout_individual_bam = bamout_individual_dir+chrom+"."+start+"."+end+".bam"
    bamout_individual_bai = bamout_individual_dir+chrom+"."+start+"."+end+".bai"
    completed = compl_folder+"SplitBamOutRegion_"+individual+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

    gwf.target_from_template("SplitBamOutRegion_"+individual+"_"+chrom+"_"+start+"_"+end,
            split_bamout_region(read_group = read_group, bamout_individual_bam = bamout_individual_bam, bamout_individual_bai = bamout_individual_bai, completed = completed))

for individual in set(itertools.chain(*GORILLA_FAMILIES)):
  for chrom, start, end in GORILLA_REGIONS:

    specie = "gorillas"
    bamout_dir = new_gatk_files+"bamout/"+specie+"/"
    bamout = bamout_dir+chrom+"."+start+"."+end+".bam"
    read_group = read_group_dir+individual+".rg.txt"

    bamout_individual_dir = bamout_files_dir+individual+"/"

    if not os.path.exists(bamout_individual_dir):
      os.makedirs(bamout_individual_dir)

    bamout_individual_bam = bamout_individual_dir+chrom+"."+start+"."+end+".bam"
    bamout_individual_bai = bamout_individual_dir+chrom+"."+start+"."+end+".bai"
    completed = compl_folder+"SplitBamOutRegion_"+individual+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

    gwf.target_from_template("SplitBamOutRegion_"+individual+"_"+chrom+"_"+start+"_"+end,
            split_bamout_region(read_group = read_group, bamout_individual_bam = bamout_individual_bam, bamout_individual_bai = bamout_individual_bai, completed = completed))

for individual in set(itertools.chain(*ORANGUTAN_FAMILIES)):
  for chrom, start, end in ORANGUTAN_AUTOSOME_REGIONS:

    specie = "orangutans"
    bamout_dir = new_gatk_files+"bamout/"+specie+"/"
    bamout = bamout_dir+chrom+"."+start+"."+end+".bam"
    read_group = read_group_dir+individual+".rg.txt"

    bamout_individual_dir = bamout_files_dir+individual+"/"

    if not os.path.exists(bamout_individual_dir):
      os.makedirs(bamout_individual_dir)

    bamout_individual_bam = bamout_individual_dir+chrom+"."+start+"."+end+".bam"
    bamout_individual_bai = bamout_individual_dir+chrom+"."+start+"."+end+".bai"
    completed = compl_folder+"SplitBamOutRegion_"+individual+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

    gwf.target_from_template("SplitBamOutRegion_"+individual+"_"+chrom+"_"+start+"_"+end,
            split_bamout_region(read_group = read_group, bamout_individual_bam = bamout_individual_bam, bamout_individual_bai = bamout_individual_bai, completed = completed))

for child, father, mother in ORANGUTAN_FAMILIES:
    for chrom in ORANGUTAN_CHROMOSOMES:

      individual = child

      bamout_individual_dir = bamout_files_dir+individual+"/"
      bam = bamout_individual_dir+chrom+".bam"
      bai = bamout_individual_dir+chrom+".bam.bai"
      completed = compl_folder+"GetIndexBam_"+individual+"_"+chrom+".COMPLETED"

      gwf.target_from_template("GetIndexBamChroms_"+individual+"_"+chrom,
              samtools_index_bam(bam, bai, completed))

for child, father, mother in GORILLA_FAMILIES:
    for chrom, start, end in GORILLA_REGIONS:

      individual = child

      bamout_individual_dir = bamout_files_dir+individual+"/"
      bam = bamout_individual_dir+chrom+"."+start+"."+end+".bam"
      bai = bamout_individual_dir+chrom+"."+start+"."+end+".bam.bai"
      completed = compl_folder+"GetIndexBam_"+individual+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

      gwf.target_from_template("GetIndexBamRegions_"+individual+"_"+chrom+"_"+start+"_"+end,
              samtools_index_bam(bam, bai, completed))

for child, father, mother in CHIMP_FAMILIES:
    for chrom, start, end in CHIMP_REGIONS:

      individual = child

      bamout_individual_dir = bamout_files_dir+individual+"/"
      bam = bamout_individual_dir+chrom+"."+start+"."+end+".bam"
      bai = bamout_individual_dir+chrom+"."+start+"."+end+".bam.bai"
      completed = compl_folder+"GetIndexBam_"+individual+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

      gwf.target_from_template("GetIndexBamRegions_"+individual+"_"+chrom+"_"+start+"_"+end,
              samtools_index_bam(bam, bai, completed))

      filtered_bam_father = bamout_files_dir+father+chrom+"."+start+"."+end+".bam"
      filtered_bam_mother = bamout_files_dir+mother+chrom+"."+start+"."+end+".bam"
      filtered_bam_child = bamout_files_dir+child+chrom+"."+start+"."+end+".bam"

      family_coverage_child_dir = family_coverage_dir+child+"/"
      if not os.path.exists(family_coverage_child_dir):
        os.makedirs(family_coverage_child_dir)

      family_coverage = family_coverage_child_dir+chrom+"."+start+"."+end+".txt"
      completed = compl_folder+"GetFamilyCoverageRegion_"+father+"_"+mother+"_"+child+"_"+chrom+"_"+start+"_"+end+".COMPLETED"

      gwf.target_from_template("GetFamilyCoverageRegion_"+father+"_"+mother+"_"+child+"_"+chrom+"_"+start+"_"+end,
                get_family_coverage_region(filtered_bam_father, filtered_bam_mother, filtered_bam_child, family_coverage, completed))

      #?? target('combine_coverage_region_' + child + '_ autosomes') << \
          #combine_coverage_region(child=child, regions=CHIMP_AUTOSOME_REGIONS, outname='autosomes')

minQ =20
minq=10

for child, father, mother in CHIMP_FAMILIES:

  chromosomes = CHIMP_AUTOSOMES
  outname = "autosomes"

  input_list = ""
  new_family_coverage = family_coverage_dir+child+"_minQ"+minQ+"_minq"+minq+"_type_"+outname".txt"
  completed = compl_folder+"CombineCoverageContextWr_"+child+"_minQ"+minQ+"_minq"+minq+"_type_"+outname+".COMPLETED"

  gwf.target_from_template("CombineCoverageContextWr_"+child+"_minQ"+minQ+"_minq"+minq+"_type_"+outname,
          call_combine_coverage_context_wr(input_list=input_list, child=child, minQ=minQ, minq=minq, chromosomes=chromosomes, new_family_coverage=new_family_coverage, completed=completed))

  for chrom in CHIMP_AUTOSOMES:

    bamout_father_dir = bamout_files_dir+father+"/"
    filtered_bam_father = bamout_father_dir+chrom+"."+start+"."+end+"."+".bam"
    filtered_bai_father = bamout_father_dir+chrom+"."+start+"."+end+"."+".bam.bai"
    filtered_mother_dir = bamout_files_dir+mother+"/"
    filtered_bam_mother = bamout_mother_dir+chrom+"."+start+"."+end+"."+".bam"
    filtered_bai_mother = bamout_mother_dir+chrom+"."+start+"."+end+"."+".bam.bai"
    filtered_child_dir = bamout_files_dir+child+"/"
    filtered_bam_child = bamout_child_dir+chrom+"."+start+"."+end+"."+".bam"
    filtered_bam_child = bamout_child_dir+chrom+"."+start+"."+end+"."+".bam.bai"
    family_coverage_child_dir = family_coverage_dir+child+"/"
    family_coverage = family_coverage_dir+chrom+"."+start+"."+end+".txt"
    completed = compl_folder+"GetFamilyCoverageContextWr_"+child+"_"+chrom+"_"+start+"_"+end+".bam"

        for chrom in CHIMP_AUTOSOMES:
        target('family_depth_chimp_' + chrom +'_' + child) << \
            get_family_coverage_context_wr(father=father, mother=mother, child=child,
                                           minQ=minQ, minq=10, refname=CHIMP_REF, chrom=chrom)


get_family_coverage_context_wr(filtered_bam_father, filtered_bam_mother, filtered_bam_child, family_coverage, completed)

        target('family_depth_chimp_' + chrom +'_' + child) << \
            get_family_coverage_context_wr(father=father, mother=mother, child=child,
                                           minQ=minQ, minq=10, refname=CHIMP_REF, chrom=chrom)
filtered_bam_father, filtered_bai_father, filtered_bam_mother, filtered_bai_mother, filtered_bam_child, filtered_bai_child, family_coverage, minQ, minq, refname, completed

def get_family_coverage_context_wr(filtered_bam_father, filtered_bai_father, filtered_bam_mother, filtered_bai_mother, filtered_bam_child, filtered_bai_child, family_coverage, minQ, minq, refname, completed)

for child, father, mother in GORILLA_FAMILIES:
    target('combine_coverage_' + child + '_ autosomes') << \
        combine_coverage_context_wr(child=child, chromosomes=GORILLA_AUTOSOMES,
                                    outname='autosomes', minQ=minQ, minq=10)
    for chrom in GORILLA_AUTOSOMES:
        target('family_depth_gorilla_' + chrom +'_' + child) << \
            get_family_coverage_context_wr(father=father, mother=mother, child=child,
                                           minQ=minQ, minq=10, refname=GORILLA_REF, chrom=chrom)

for child, father, mother in ORANGUTAN_FAMILIES:
    target('combine_coverage_' + child + '_ autosomes') << \
        combine_coverage_context_wr(child=child, chromosomes=ORANGUTAN_AUTOSOMES,
                                    outname='autosomes', minQ=minQ, minq=10)
    for chrom in ORANGUTAN_AUTOSOMES:
        target('family_depth_orangutan_' + chrom +'_' + child) << \
            get_family_coverage_context_wr(father=father, mother=mother, child=child,
                                           minQ=minQ, minq=10, refname=ORANGUTAN_REF, chrom=chrom)

get_family_coverage_region = \
    template(input=['bamout_files/{father}/{chrom}.{start}.{end}.bam',
                    'bamout_files/{mother}/{chrom}.{start}.{end}.bam',
                    'bamout_files/{child}/{chrom}.{start}.{end}.bam'],
             output=['new_family_coverage/{child}/{chrom}.{start}.{end}.txt'],
             walltime='10:00:00',
             memory='1g',
             account="DanishPanGenome") << '''
source /com/extra/samtools/LATEST/load.sh
mkdir -p new_family_coverage/{child}/
samtools depth bamout_files/{father}/{chrom}.{start}.{end}.bam bamout_files/{mother}/{chrom}.{start}.{end}.bam bamout_files/{child}/{chrom}.{start}.{end}.bam | awk '$3>=5 && $3<=150 && $4>=5 && $4<=150 && $5>=5 && $5<=150' | cut -f3- | ~/Scripts/table.py > new_family_coverage/{child}/{chrom}.{start}.{end}.txt
'''

