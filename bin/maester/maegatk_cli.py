#!/usr/bin/python

import os
import click
import os.path
import sys
import pysam
import pandas as pd
from maegatk.maegatk.maegatkHelp import *
from ruamel import yaml
from ruamel.yaml.scalarstring import SingleQuotedScalarString as sqs
from multiprocessing import Pool


##


def main():


    # Args

    # Optional args. STRICTLY, THE ONLY ONE THAT NEED TO BE PASSED.
    script_dir = sys.argv[1]
    bam = sys.argv[2]
    ncores = sys.argv[3]
    barcodes = sys.argv[4]
    min_reads = sys.argv[5]
    min_base_qual = sys.argv[6]
    min_alignment_quality = sys.argv[7]
    output = os.getcwd()

    # Default args
    mito_genome = 'rCRS'
    barcode_tag = 'CR'
    min_barcode_reads = 100
    nsamples = 1500
    umi_barcode = 'UR'
    nhmax = 2
    nmmax = 15
    max_javamem = '6000m'
    skip_r = True
    jobs = 0
    name = 'maegatk'
    snake_stdout = True
    cluster = ''


    ##


    # CHECK-IN: CORES, REFERENCE GENOME, BAM, BARCODES
    
    ## Cores
    if ncores == "detect":
        ncores = str(available_cpu_count())
    else:
        ncores = str(ncores)
    
        
    ##
    
    
    ## .bam
    filename, file_extension = os.path.splitext(bam)
    
    if(file_extension != ".bam"):
        sys.exit('ERROR: in `bcall` mode, the input should be an individual .bam file.')
    if not os.path.exists(bam):
        sys.exit('ERROR: No file found called "' + bam + '"; please specify a valid .bam file.')
    if not os.path.exists(bam + ".bai"):
        sys.exit('ERROR: index your input .bam file for `bcall` mode.')
        
    if barcode_tag == "X":
        sys.exit('ERROR: in `bcall` mode, must specify a valid read tag ID (generally two letters).')
    
    click.echo(gettime() + "Found bam file: " + bam + " for genotyping.")	
    
    
    ##
    
    
    ## Reference
    of = output
    tf = of + "/temp" 
    bcbd = tf + "/barcoded_bams" 
    folders = [ of, tf, bcbd, of + "/final" ] 
        
    for x in folders:
        make_folder(x) 
        
    rawsg = os.popen('ls ' + script_dir + '/bin/anno/fasta/*.fasta').read().strip().split("\n") 
    supported_genomes = [ 
        x.replace(script_dir + "/bin/anno/fasta/", "").replace(".fasta", "") \
        for x in rawsg 
    ]  
    
    fastaf, mito_chr, mito_length = handle_fasta_inference(mito_genome, supported_genomes, script_dir, 'bcall', of)
    idxs = pysam.idxstats(bam).split("\n") 
    
    try:
        bam_length = [ int(x.split('\t')[1]) for x in idxs if x.startswith(mito_chr) ][0]
    except:
        bam_length = 0
    
    if(mito_length == bam_length):
        click.echo(gettime() + "User specified mitochondrial genome matches .bam file")
    elif(bam_length == 16569):
        click.echo(gettime() + "User specified mitochondrial genome does NOT match .bam file; using rCRS instead (length == 16569)")
        fastaf, mito_chr, mito_length = handle_fasta_inference("rCRS", supported_genomes, script_dir, 'bcall', of)
    elif(bam_length == 16571):
        click.echo(gettime() + "User specified mitochondrial genome does NOT match .bam file; using hg19 instead (length == 16571)")
        fastaf, mito_chr, mito_length = handle_fasta_inference("hg19", supported_genomes, script_dir, 'bcall', of)
    else:
        click.echo(gettime() + "User specified mitochondrial genome does NOT match .bam file; correctly specify reference genome or .fasta file")
        quit()
    
        
    ##
    
    
    ## Barcodes 
    if (os.path.exists(barcodes)): 
        find_barcodes = False
    else:
        find_barcodes = True
        
    if find_barcodes:
        barc_quant_file = of + "/final/barcodeQuants.tsv"
        passing_barcode_file = of + "/final/passingBarcodes.tsv"
        find_barcodes_py = script_dir + "/bin/python/find_barcodes.py"
        
        pycall = " ".join(['python', find_barcodes_py, bam, bcbd, barcode_tag, str(min_barcode_reads), mito_chr, barc_quant_file, passing_barcode_file])
        os.system(pycall)
        barcodes = passing_barcode_file
    
    barcode_files = split_barcodes_file(barcodes, nsamples, output)
    split_barcoded_bam_py = script_dir + "/bin/python/split_barcoded_bam.py"
    
    for i in range(len(barcode_files)):
        one_barcode_file = barcode_files[i]
        pycall = " ".join(['python', split_barcoded_bam_py, bam, bcbd, barcode_tag, one_barcode_file, mito_chr])
        os.system(pycall)
        
    click.echo(gettime() + "Finished determining/splitting barcodes for genotyping.")
        
    
    ##
    
    
    ## Summary input checks
    click.echo(gettime() + f'Input .bam: {bam}')
    click.echo(gettime() + f'Reference mito genome: {mito_genome}')
    click.echo(gettime() + f'Min reads: {min_reads}')
    click.echo(gettime() + f'Find barcodes: {find_barcodes}')
    click.echo(gettime() + f'N cores: {ncpus}')
    
    
    ##
    
    
    # CHECK-IN pt II: discard low quality cells
    bam = bcbd
    
    if len(os.listdir(bam)) == 0:
        sys.exit('ERROR: Could not import any samples from the user specification; check flags, logs and input configuration; QUITTING')
    else:
        bams = [ f'{bam}/{x}' for x in os.listdir(bcbd) ]
    
    samples = []
    samplebams = []
        
    for bam in bams:
        base = os.path.basename(bam)
        basename = os.path.splitext(base)[0]
        samples.append(basename)
        samplebams.append(bam)
    
    pool = Pool(processes=int(ncores))
    pm = pool.map(verify_bai, samplebams)
    pool.close()
    
    samples_fail = []
    for i in range(len(samples)):
        sample = samples[i]
        bam = samplebams[i]
        if not verify_sample_mitobam(bam, mito_chr, mito_length):
            samples_fail.append(sample)  
    
    if len(samples_fail) > 0:
        click.echo(gettime() + "NOTE: There are failed samples...")
        rmidx = findIdx(samples, samples_fail)
        for index in sorted(rmidx, reverse=True):
            print("REMOVED: ", samples[index])
            del samples[index]
            del samplebams[index]
            
            
    ## 
    
    
    ## Exit if None remaining...      
    if not len(samples) > 0:
        sys.exit('ERROR: Could not import any samples from the user specification. \nERROR: check flags, logs, and input configuration (including reference mitochondrial genome); \nQUITTING')
    
    
    ##
    
    
    # PREP SNAKEMAKE FOLDERS AND log
    
    # Folders
    of = output
    tf = of + "/temp"
    qc = of + "/qc"
    logs = of + "/logs"
    
    folders = [ 
        logs, 
        of + "/logs/filterlogs", 
        of + "/logs/rmdupslogs",
        of + "/fasta", 
        of + "/.internal",
        of + "/.internal/parseltongue",
        of + "/.internal/samples", of + "/final", 
        tf, 
        tf + "/ready_bam", 
        tf + "/temp_bam",
        tf + "/sparse_matrices", 
        tf + "/quality",
        qc, 
        qc + "/quality", 
        qc + "/depth"
    ]
    
    for x in folders:
        make_folder(x)
    
        
    ##
    
    
    # logs
    logf = open(output + "/logs" + "/base.maegatk.log", 'a')
    
    # Internals...
    if not os.path.exists(of + "/.internal/README"):
        with open(of + "/.internal/README" , 'w') as outfile:
            outfile.write("This folder creates important (small) intermediate; don't modify it.\n\n")
    if not os.path.exists(of + "/.internal/parseltongue/README"):	
        with open(of + "/.internal/parseltongue/README" , 'w') as outfile:
            outfile.write("This folder creates intermediate output to be interpreted by Snakemake; don't modify it.\n\n")
    if not os.path.exists(of + "/.internal/samples/README"):
        with open(of + "/.internal" + "/samples" + "/README" , 'w') as outfile:
            outfile.write("This folder creates samples to be interpreted by Snakemake; don't modify it.\n\n")
    
    for i in range(len(samples)):
        with open(of + "/.internal/samples/" + samples[i] + ".bam.txt" , 'w') as outfile:
            outfile.write(samplebams[i])
    
    
    ##
    
    
    # HERE WE GO: SCATTER!
    
    click.echo(gettime() + f"Starting analysis with maegatk: {len(samples)} to process...", logf)
    
    # Conf and call
    dict1 = {
        'input_directory' : sqs(input), 'output_directory' : sqs(output), 'script_dir' : sqs(script_dir),
        'fasta_file' : sqs(fastaf), 'mito_chr' : sqs(mito_chr), 'mito_length' : sqs(mito_length), 
        'base_qual' : sqs(min_base_qual), 'umi_barcode' : sqs(umi_barcode),
        'alignment_quality' : sqs(min_alignment_quality), 
        'NHmax' : sqs(nhmax), 'NMmax' : sqs(nmmax), 'min_reads' : sqs(min_reads),'max_javamem' : sqs(max_javamem)
    }
    
    snakeclust = ""
    njobs = int(jobs)
    
    if njobs > 0 and cluster != "":
        snakeclust = " --jobs " + jobs + " --cluster '" + cluster + "' " 
        click.echo(gettime() + "Recognized flags to process jobs on a computing cluster.", logf)
        click.echo(gettime() + "Processing samples with "+ ncores +" threads", logf) 
    
    y_s = of + "/.internal/parseltongue/snake.scatter.yaml"
    with open(y_s, 'w') as yaml_file:
        yaml.dump(dict1, yaml_file, default_flow_style=False, Dumper=yaml.RoundTripDumper)
    
    
    snake_stats = logs + "/" + name + ".snakemake_scatter.stats"
    snake_log = logs + "/" + name + ".snakemake_scatter.log"
    
    snake_log_out = ""
    
    if not snake_stdout:
        snake_log_out = ' &>' + snake_log
    
    snakecmd_scatter = 'snakemake' + snakeclust + ' --snakefile ' + script_dir + '/bin/snake/Snakefile.maegatk.Scatter --cores '+ ncores +' --config cfp="'  + y_s + '" --stats '+snake_stats + snake_log_out 
    click.echo(gettime() + "OK until scatter!!")
    
    # Run 
    os.system(snakecmd_scatter)
    click.echo(gettime() + "OK post scatter!!")
    
    
    ##
    
    
    # HERE WE GO: GATHER!
    
    maegatk_directory = output
    
    dict2 = {
        'maegatk_directory' : sqs(maegatk_directory), 'name' : sqs(name), 'script_dir' : sqs(script_dir)
    }
    
    y_g = maegatk_directory + "/.internal/parseltongue/snake.gather.yaml"
    with open(y_g, 'w') as yaml_file:
        yaml.dump(dict2, yaml_file, default_flow_style=False, Dumper=yaml.RoundTripDumper)
    
    snake_stats = logs + "/" + name + ".snakemake_gather.stats"
    snake_log = logs + "/" + name + ".snakemake_gather.log"
    
    snakecmd_gather = 'snakemake --snakefile ' + script_dir + '/bin/snake/Snakefile.maegatk.Gather --cores '+ncores+' --config cfp="' + y_g + '" --stats '+snake_stats + snake_log_out 
    click.echo(gettime() + "OK until gather!!")
    
    # Run
    os.system(snakecmd_gather)
    click.echo(gettime() + "OK post gather!!")
        click.echo(gettime() + "Successfully created final output files", logf)
    
    
    ##
    
    
# Run
if __name__ == '__main__':
    main()


##
