from gwf import Workflow

gwf = Workflow()

refGenome = 'panTro5'

#Directory containing reference genome'
refdir = 'ref-genomes'

#index fasta file
gwf.target(
  f'index_fastafile',
  inputs=[f'{refdir}/{refGenome}.fa'],
  outputs=[f'{refdir}/{refGenome}.amb', f'{refdir}/{refGenome}.ann',
           f'{refdir}/{refGenome}.pac', f'{refdir}/{refGenome}.sa'],
  walltime="12:00:00", memory="8g") << \
  f'''
  bwa index -p {refdir}/{refGenome} -a bwtsw {refdir}/{refGenome}.fa
  '''

f = open('fastq_files.txt')
name2fastq_pairs = {}
for line in f:
    name, fq1, fq2 = line.split()
    if name.lower() == 'individual':
        continue
    if name not in name2fastq_pairs:
        name2fastq_pairs[name] = []
    name2fastq_pairs[name].append((fq1, fq2))

for name in name2fastq_pairs:
    bam_files = []
    for i in range(len(name2fastq_pairs[name])):
        R1, R2 = name2fastq_pairs[name][i]
        bamfile = f'bam_files/{name}_{i}.bam'
        bam_files.append(bamfile)

        #Map reads and produce sorted bam files with duplicates removed.
        gwf.target(
            f'bwa_map_{name}_{i}',
            inputs=[R1, R2, f'{refdir}/{refGenome}.amb', f'{refdir}/{refGenome}.ann',
                    f'{refdir}/{refGenome}.pac', f'{refdir}/{refGenome}.sa'],
            outputs=[bamfile],
            cores=16, walltime="100:00:00", memory="64g") << \
            f'''
            bwa mem -t 16 {refdir}/{refGenome} {R1} {R2} | \
            samtools sort | \
            samtools rmdup -s - {bamfile}
            '''

    inputbams = ' '.join(bam_files)
    
    #Merge bamfiles
    gwf.target(
        f'merge_bam_files_{name}',
        inputs=bam_files,
        outputs=[f'merged_bam_files/{name}.rg.txt', f'merged_bam_files/{name}.bam'],
        walltime="100:00:00", memory="8g") << \
        f'''
        rg_file=merged_bam_files/{name}.rg.txt

        rm -f $rg_file
        for input in {inputbams}; do
            bn=`basename $input .bam`
            echo -e "@RG\tID:$bn\tSM:{name}\tLB:{name}\tPL:Illumina" >> $rg_file
        done

        samtools merge -rh $rg_file - {inputbams} | samtools rmdup -s - merged_bam_files/{name}.bam
        '''
