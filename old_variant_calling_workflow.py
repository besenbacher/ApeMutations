from gwf import *
from gwf_templates import *
import itertools

### Workflow for variant calling. ###

ORANGUTAN_REF = 'ponAbe2'
CHIMP_REF     = 'panTro5' 
GORILLA_REF   = 'gorGor4'

ORANGUTANS = ['Buschi', 'Moni', 'Masala', 'Farida', 'Ogan',  'Aris', 'Schubbi', 'Pongo']
CHIMPS     = ['Frits', 'Carolina', 'Carl', 'Simliki', 'ERR466113','ERR466114','ERR466115','ERR466116','ERR466117','ERR466118','ERR466119','ERR466120','ERR466121']
GORILLAS   = ['Banjo', 'Mimi', 'Mutasi', 'Mawenzi', 'Efata', 'Undi']

def read_family_description(fname):
    f = open(fname)
    res = []
    for line in f:
        L = line.split()
        res.append((L[0],L[1],L[2]))
    f.close()
    return res

CHIMP_FAMILIES = read_family_description('family_description/chimpanzees.txt')
GORILLA_FAMILIES = read_family_description('family_description/gorillas.txt')
ORANGUTAN_FAMILIES = read_family_description('family_description/orangutans.txt')

#Order: (Child, Father, Mother)

ALL_APES = ORANGUTANS + CHIMPS + GORILLAS

#First we need to index the reference genomes with Picard
#refGenomes should be in ./ref-genomes

target('picard_index_orang')   << picard_index_reference(refGenome=ORANGUTAN_REF)
target('picard_index_chimp')   << picard_index_reference(refGenome=CHIMP_REF)
target('picard_index_gorilla') << picard_index_reference(refGenome=GORILLA_REF)

## Use GATK to process bam files before calling.

#To work with the bam files in GATK we first need them indexed...
# bamfiles should be in ./merged-bams

for individual in ALL_APES:
    target('samtools_index_bam_' + individual) << samtools_index_bam(individual=individual)
    target('samtools_index_filtered_bam_' + individual) << samtools_index_filtered_bam(individual=individual)


#Next, identify regions that needs to be locally re-aligned.
# intermediary output files from gatk will be put in ./gatk_files/

for individual in CHIMPS:
    target('collect_realign_regions_' + individual) << collect_realign_regions(refGenome=CHIMP_REF, individual=individual)

for individual in ORANGUTANS:
    target('collect_realign_regions_' + individual) << collect_realign_regions(refGenome=ORANGUTAN_REF, individual=individual)

for individual in GORILLAS:
    target('collect_realign_regions_' + individual) << collect_realign_regions(refGenome=GORILLA_REF, individual=individual)


#Realign the bam files:

for individual in CHIMPS:
    target('realign_' + individual) << realign(refGenome=CHIMP_REF, individual=individual)

for individual in ORANGUTANS:
    target('realign_' + individual) << realign(refGenome=ORANGUTAN_REF, individual=individual)

for individual in GORILLAS:
    target('realign_' + individual) << realign(refGenome=GORILLA_REF, individual=individual)


#Recalibrate base quality scores:

orang_known = ['../KnownVariants/abelii.2012-08-17.snps.dropsamples.autos.vqsr99.vcf',
               '../KnownVariants/abelii.2012-08-17.snps.dropsamples.chrX.vcf',
               '../KnownVariants/pygmaeus.2012-08-17.snps.dropsamples.autos.vqsr99.vcf',
               '../KnownVariants/pygmaeus.2012-08-17.snps.dropsamples.chrX.vcf']

chimp_known = ['../faststorage/KnownVariants/pantro.snps.2012-11-25.autos.99vqsr.beagle_' + CHIMP_REF +'_lift.vcf',
               '../faststorage/KnownVariants/gatk.allchimpanzee.2012-05-30.snps.chrX.vqsr99_' + CHIMP_REF + '_lift.vcf']
gorilla_known = ['../faststorage/KnownVariants/emitall-var.dropsamples.autos.vqsr99.2012-04-26.abfilter_' + GORILLA_REF  + '_lift_sorted.vcf']

for individual in CHIMPS:
    target('calc_recalibrate_info_' + individual) << \
        calc_recalibrate_info(individual=individual,
                              refGenome=CHIMP_REF,
                              knownlist=chimp_known)
    target('recalibrate_' + individual) << recalibrate(refGenome=CHIMP_REF, individual=individual)
    target('filter_bam_' + individual) << filter_bam_files(individual=individual,
                                                           refGenome=CHIMP_REF)


for individual in ORANGUTANS:
    target('calc_recalibrate_info_' + individual) << \
        calc_recalibrate_info(individual=individual,
                              refGenome=ORANGUTAN_REF,
                              knownlist=orang_known)
    target('recalibrate_' + individual) << recalibrate(refGenome=ORANGUTAN_REF, individual=individual)
    target('filter_bam_' + individual) << filter_bam_files(individual=individual,
                                                           refGenome=ORANGUTAN_REF)

for individual in GORILLAS:
    target('calc_recalibrate_info_' + individual) << \
        calc_recalibrate_info(individual=individual,
                              refGenome=GORILLA_REF,
                              knownlist=gorilla_known)
    target('recalibrate_' + individual) << recalibrate(refGenome=GORILLA_REF, individual=individual)
    target('filter_bam_' + individual) << filter_bam_files(individual=individual,
                                                          refGenome=GORILLA_REF)


#Call variants using haplotype caller

def get_chromosomes(refname):
    f = open('ref-genomes/'+ refname + '.fa.fai')
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
    f = open('split_files/' + refname + '_split_1000_1000000.txt')
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
  target('hapcaller_bo_orangutan_' + chrom + '_' + str(start)) << \
    haplotype_caller_region_bo(refGenome=ORANGUTAN_REF, specie="orangutans",
                               individuals=ORANGUTANS, chrom=chrom, start=start, end=end)


for chrom, start, end in GORILLA_REGIONS:
    target('hapcaller_bo_gorilla_' + chrom + '_' + str(start)) << \
      haplotype_caller_region_bo(refGenome=GORILLA_REF, specie="gorillas",
                                individuals=GORILLAS, chrom=chrom, start=start, end=end)


for chrom, start, end in CHIMP_REGIONS:
   target('hapcaller_bo_chimps_' + chrom + '_' + str(start)) << \
     haplotype_caller_region_bo(refGenome=CHIMP_REF, specie="chimpanzees",
                               individuals=CHIMPS, chrom=chrom, start=start, end=end)


target('combine_regions_chimps') << combine_regions(refname=CHIMP_REF, specie="chimpanzees", regions=CHIMP_REGIONS)
target('combine_regions_gorillas') << combine_regions(refname=GORILLA_REF, specie="gorillas", regions=GORILLA_REGIONS)
target('combine_regions_orangutans') << combine_regions(refname=ORANGUTAN_REF, specie="orangutans", regions=ORANGUTAN_AUTOSOME_REGIONS)

for individual in ALL_APES:
    target('get_rg_' + individual) << get_readgroup(individual=individual)

for individual in set(itertools.chain(*CHIMP_FAMILIES)):
    for chrom, start, end in CHIMP_REGIONS:
        target('split_bo_' +individual + '_' +chrom +'_' + str(start)) << \
          split_bamout_region(individual=individual, chrom=chrom, start=start, end=end, specie="chimpanzees")

for individual in set(itertools.chain(*GORILLA_FAMILIES)):
  for chrom, start, end in GORILLA_REGIONS:
      target('split_bo_' +individual + '_' +chrom +'_' + str(start)) << \
          split_bamout_region(individual=individual, chrom=chrom, start=start, end=end, specie="gorillas")

for individual in set(itertools.chain(*ORANGUTAN_FAMILIES)):
    for chrom, start, end in ORANGUTAN_AUTOSOME_REGIONS:
        target('split_bo_' +individual + '_' +chrom +'_' + str(start)) << \
            split_bamout_region(individual=individual, chrom=chrom, start=start, end=end, specie="orangutans")

for child, father, mother in ORANGUTAN_FAMILIES:
    for chrom in ORANGUTAN_CHROMOSOMES:
        target('samtools_index_bamout_' + child + '_' + chrom) << \
            samtools_index_bamout_chrom(individual=child, chrom=chrom)

for child, father, mother in GORILLA_FAMILIES:
    for chrom, start, end in GORILLA_REGIONS:
        target('samtools_index_bamout_' + child + '_' + chrom +'_' + str(start) ) << \
            samtools_index_bamout(individual=child, chrom=chrom, start=start, end=end)

for child, father, mother in CHIMP_FAMILIES:
    for chrom, start, end in CHIMP_REGIONS:
        target('samtools_index_bamout_' + child + '_' + chrom +'_' + str(start) ) << \
            samtools_index_bamout(individual=child, chrom=chrom, start=start, end=end)
       target('get_coverage_region_' + child + '_' + chrom +'_' + str(start)) << \
           get_family_coverage_region(child=child, father=father, mother=mother, chrom=chrom, start=start, end=end, outname='autosomes')
   target('combine_coverage_region_' + child + '_ autosomes') << \
       combine_coverage_region(child=child, regions=CHIMP_AUTOSOME_REGIONS, outname='autosomes')


minQ =20
for child, father, mother in CHIMP_FAMILIES:
    target('combine_coverage_' + child + '_ autosomes') << \
        combine_coverage_context_wr(child=child, chromosomes=CHIMP_AUTOSOMES,
                                    outname='autosomes', minQ=minQ, minq=10)
    for chrom in CHIMP_AUTOSOMES:
        target('family_depth_chimp_' + chrom +'_' + child) << \
            get_family_coverage_context_wr(father=father, mother=mother, child=child,
                                           minQ=minQ, minq=10, refname=CHIMP_REF, chrom=chrom)

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