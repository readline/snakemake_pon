import os
import json

snakedir = os.getcwd()
configfile: 'config.yaml'

print('Snakemake running dir:', snakedir)
print('Pipeline working dir:', config['workdir'])
print('#'*100,'\n')

print('Config setted to:')
print(json.dumps(config, indent=4, sort_keys=True), '\n')
print('#'*100,'\n')

bamdic = ['%s/%s/%s'%(config['workdir'],config['bamsubfolder'],i) for i in os.listdir('%s/%s'%(config['workdir'],config['bamsubfolder'])) if i[-4:] in ['.bam','cram']]
bamdic = {i.split('/')[-1].split('.')[0]:i for i in bamdic}
print('No. of samples used for this PoN:',len(bamdic))
for i in bamdic:
    print('>>> %s: %s'%(i,bamdic[i]))

itvlist = ['itv_%.4d'%i for i in range(1,config['references']['interval']+1)]
print('Distribute to intervals:',len(itvlist))

workdir: config['workdir']

rule all:
    input:
        mutect = 'Mutect2.pon.vcf.gz',
        cnvkit = 'CNVkit.pon.cnn',

rule Prepare:
    input:
    output:
        itv = snakedir+'/ref/intervals/itv.17.ok',
    log:
        out = snakedir+'/logs/A1.Prepare/All.o',
        err = snakedir+'/logs/A1.Prepare/All.e',
    threads:  2
    resources:
        mem = '4g',
        extra = ' --gres=lscratch:10 ',
    run:
        shell(
            'rm -rf {snakedir}/singularity\n'
            'mkdir -p {snakedir}/singularity\n'
            'cd {snakedir}/ref\n'
            'rm -f *.fa *.fa.fai *.dict *.bed af-only-gnomad.hg38.vcf.gz*'
            ' > {log.out} 2> {log.err}\n'
            'ln -s {config[references][fasta]} ref.fa'
            ' >> {log.out} 2>> {log.err}\n'
            'ln -s {config[references][fastafai]} ref.fa.fai'
            ' >> {log.out} 2>> {log.err}\n'
            'ln -s {config[references][fastadict]} ref.dict'
            ' >> {log.out} 2>> {log.err}\n'
            'wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz\n'
            'wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi\n'
            'wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz \n'
            'gunzip refFlat.txt.gz \n'
             )

        # Prepare singularity containers
        for simgname in config['simg']:
            shell(
                'cd {snakedir}/singularity\n'+\
                'wget {config[simg][%s]} -O %s.simg'%(simgname, simgname)+\
                ' >> {log.out} 2>> {log.err} &&\n'+\
                'echo {simgname}.simg: {config[simg][%s]} >> {snakedir}/singularity/simg.log'%(simgname)
            )
        # singularity container checkpoint
        if os.path.exists('{snakedir}/singularity/simg.log'):
            with open('{snakedir}/singularity/simg.log') as infile:
                print(len(infile.readlines()),infile.readlines())
                if len(infile.readlines()) != len(config['simg']):
                    shell(
                        'cd {snakedir}/singularity\n'
                        'rm {snakedir}/singularity/simg.log \n'
                        'echo Error! Singularity images numbers doesn\'t match.'
                    )
                    
        # Prepare WXS bed files
        if 'flankbed' in config['references']:
            shell(
                'cd {snakedir}/ref\n'
                'ln -s {config[references][wxsbed]} wxs.bed'
                ' >> {log.out} 2>> {log.err}\n'
                'ln -s {config[references][flankbed]} wxs.flank.bed'
                ' >> {log.out} 2>> {log.err}\n'
            )
            
        # Prepare exclusion file    
        shell(
            'cd {snakedir}/ref\n'
            'wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz'
            ' >> {log.out} 2>> {log.err}\n'
            'wget https://github.com/lh3/varcmp/raw/master/scripts/LCR-hs38.bed.gz'
            ' >> {log.out} 2>> {log.err}\n'
            'wget https://www.encodeproject.org/files/ENCFF220FIN/@@download/ENCFF220FIN.bed.gz'
            ' >> {log.out} 2>> {log.err}\n'
            'module load singularity'
            ' >> {log.out} 2>> {log.err}\n'
            'zcat ENCFF220FIN.bed.gz hg38-blacklist.v2.bed.gz LCR-hs38.bed.gz|cut -f1-3 > merge.1.bed'
            ' 2>> {log.err}\n'
            '{config[singularity]} ../singularity/bedtools.simg bedtools sort -i merge.1.bed > merge.2.bed'
            ' 2>> {log.err}\n'
            '{config[singularity]} ../singularity/bedtools.simg bedtools merge -i merge.2.bed > Merge.exclude.bed'
            ' 2>> {log.err}\n'
            'rm merge.1.bed merge.2.bed ENCFF220FIN.bed.gz hg38-blacklist.v2.bed.gz LCR-hs38.bed.gz'
        )
        
        # Unpack interval files
        if os.path.exists('{snakedir}/ref/intervals/itv.17.ok'):
            pass
        else:
            shell(
                'cd {snakedir}/ref/intervals\n'
                'tar zxvf intervals.tar.gz'
                ' >> {log.out} 2>> {log.err}\n'
            )
        

rule Mutect_call:
    input:
        idxready = snakedir+'/ref/intervals/itv.17.ok',
        bam = lambda wildcards: bamdic[wildcards.sample],
    output:
        vcf = temp('Mutect2/{sample}/{sample}.{itv}.vcf.gz'),
    log:
        out = snakedir+'/logs/B1.Mutect_call/{sample}.{itv}.o',
        err = snakedir+'/logs/B1.Mutect_call/{sample}.{itv}.e',
    threads:  4
    resources:
        mem  = '36g',
        extra = ' --gres=lscratch:50 ',
    shell:
        'module load singularity > {log.out} 2> {log.err}\n'
        '{config[singularity]} {snakedir}/singularity/gatk.simg '
        'gatk --java-options \"-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xmx32G\" '
        ' Mutect2'
        ' -R {snakedir}/ref/ref.fa'
        ' -I {input.bam}'
        ' --max-mnp-distance 0'
        ' -L {snakedir}/ref/intervals/{wildcards.itv}.bed'
        ' -O {output.vcf}'
        ' >> {log.out} 2>> {log.err}\n'
            

rule Mutect_pon:
    input:
        vcfs = lambda wildcards: expand('Mutect2/{sample}/{sample}.{itv}.vcf.gz', sample=bamdic.keys(), itv=wildcards.itv),
    output:
        gdb = directory('Mutect2/{itv}.gdb'),
        vcf = temp('Mutect2/{itv}.pon.vcf.gz'),
        vcfi = temp('Mutect2/{itv}.pon.vcf.gz.tbi'),
    log:
        out = snakedir+'/logs/B2.Mutect_pon/{itv}.o',
        err = snakedir+'/logs/B2.Mutect_pon/{itv}.e',
    threads:  4
    resources:
        mem  = '48g',
        extra = ' --gres=lscratch:50 ',
    run:
        vcfcmd = ''
        for vcf in input.vcfs:
            vcfcmd += ' -V %s'%(vcf)
        shell(
            'rm -rf {output.gdb}\n'
            'module load singularity > {log.out} 2> {log.err}\n'
            '{config[singularity]} {snakedir}/singularity/gatk.simg '
            'gatk --java-options \"-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xmx46G\" '
            ' GenomicsDBImport'
            ' -R {snakedir}/ref/ref.fa'
            ' -L {snakedir}/ref/intervals/{wildcards.itv}.bed'
            ' {vcfcmd}'
            ' --genomicsdb-workspace-path {output.gdb}'
            ' >> {log.out} 2>> {log.err}\n'
            '{config[singularity]} {snakedir}/singularity/gatk.simg '
            'gatk --java-options \"-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xmx46G\" '
            ' CreateSomaticPanelOfNormals'
            ' -R {snakedir}/ref/ref.fa'
            ' -V gendb://{output.gdb}'
            ' -O {output.vcf}'
            ' --germline-resource {snakedir}/ref/af-only-gnomad.hg38.vcf.gz'
            ' --max-germline-probability 0.5'
            ' >> {log.out} 2>> {log.err}\n'
            'ls {output.gdb}'
            ' >> {log.out} 2>> {log.err}\n'
            'ls Mutect2'
            ' >> {log.out} 2>> {log.err}\n'
        )
            
            
rule Mutect_merge:
    input:
        vcfs = lambda wildcards: expand('Mutect2/{itv}.pon.vcf.gz', itv=itvlist),
    output:
        vcf = 'Mutect2.pon.vcf.gz',
    params:
        vcf = 'Mutect2/Mutect2.pon.vcf.gz',
    log:
        out = snakedir+'/logs/B3.Mutect_merge/Merge.o',
        err = snakedir+'/logs/B3.Mutect_merge/Merge.e',
    threads:  4
    resources:
        mem  = '36g',
        extra = ' --gres=lscratch:50 ',
    run:
        vcfcmd = ''
        for vcf in input.vcfs:
            vcfcmd += ' -I %s'%(vcf)
        shell(
            'module load singularity > {log.out} 2> {log.err}\n'
            '{config[singularity]} {snakedir}/singularity/gatk.simg '
            'gatk --java-options \"-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xmx32G\" '
            '  MergeVcfs'
            '  {vcfcmd}'
            '  -O {params.vcf}'
            '  >> {log.out} 2>> {log.err}\n'
        )
        if 'flankbed' in config['references']:
            shell(
                'module load singularity > {log.out} 2> {log.err}\n'
                '{config[singularity]} {snakedir}/singularity/gatk.simg '
                'gatk --java-options \"-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xmx12G\" '
                '  SelectVariants'
                '  -V {params.vcf}'
                '  -L {config[references][flankbed]}'
                '  -O {output.vcf}'
                '  >> {log.out} 2>> {log.err}\n'
                'rm -rf Mutect2'
                '  >> {log.out} 2>> {log.err}\n'
                )
        else:
            shell(
                'mv {params.vcf} {output.vcf}'
                '  >> {log.out} 2>> {log.err}\n'
                'mv {params.vcf}.tbi {output.vcf}.tbi'
                '  >> {log.out} 2>> {log.err}\n'
                'rm -rf Mutect2'
                '  >> {log.out} 2>> {log.err}\n'
                )
        
rule Cnvkit_pon:
    input:
        idxready = snakedir+'/ref/intervals/itv.17.ok',
        
    output:
        cnn = 'CNVkit.pon.cnn',
    log:
        out = snakedir+'/logs/C1.Cnvkit_pon/prepare.o',
        err = snakedir+'/logs/C1.Cnvkit_pon/prepare.e',
    threads:  24
    resources:
        mem  = '64g',
        extra = ' --gres=lscratch:50 ',
    run:
        if 'wxsbed' in config['references']:
            wxscmd = ' -m hybrid --targets %s '%(config['references']['wxsbed'])
        bamcmd = ''
        for sample in bamdic:
            bamcmd += ' %s'%bamdic[sample]
        shell(
            'module load singularity > {log.out} 2> {log.err}\n'
            '{config[singularity]} {snakedir}/singularity/cnvkit.simg '
            '  cnvkit.py access '
            '  {snakedir}/ref/ref.fa '
            '  -x {snakedir}/ref/Merge.exclude.bed '
            '  -o {snakedir}/ref/access-excludes.hg38.bed'
            '  >> {log.out} 2>> {log.err}\n'
            '{config[singularity]} {snakedir}/singularity/cnvkit.simg '
            '  cnvkit.py batch'
            '  --fasta {snakedir}/ref/ref.fa'
            '  --access {snakedir}/ref/access-excludes.hg38.bed'
            '  --annotate {snakedir}/ref/refFlat.txt'
            '  -p {threads}'
            '  {wxscmd} '
            '  --normal {bamcmd}'
            '  --output-reference {output.cnn}'
            '  >> {log.out} 2>> {log.err}\n'
            'rm *.antitargetcoverage.cnn'
            '  >> {log.out} 2>> {log.err}\n'
            'rm *.targetcoverage.cnn'
            '  >> {log.out} 2>> {log.err}\n'
            'rm *.antitarget.bed *.target.bed'
            '  >> {log.out} 2>> {log.err}\n'
        )
