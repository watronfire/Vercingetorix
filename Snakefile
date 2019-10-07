import os
from fnmatch import fnmatch
import subprocess
import pandas as pd

# Random ass shit: mixcr and migec don't do well with spaces in file names. Maybe even parenthesis.

in_dir = config["__default__"]["in_dir"]
in_ext = config["__default__"]["in_ext"]
out_dir = config["__default__"]["out_dir"]

SAMPLES = dict()
for f in os.listdir( in_dir ):
    if fnmatch( f, "*.{}".format( in_ext ) ):
        fsplit = f.split( "_" )
        if fsplit[0] not in SAMPLES:
            SAMPLES[fsplit[0]] = [f]
        else:
            SAMPLES[fsplit[0]].append( f )

rule all:
    input:
        expand( os.path.join( out_dir, "analyze/{sample}.clonotypes.{chain}.txt" ), out_dir=out_dir, sample=SAMPLES, chain=["TRA", "TRB"] )

# TODO combine log files into a pandas dataframe and export into csv.

rule repertoire_assembly:
    input:
        os.path.join( out_dir, "cleaned/{sample}_R1.fastq.gz" ),
        os.path.join( out_dir, "cleaned/{sample}_R2.fastq.gz" )
    output:
        os.path.join( out_dir, "analyze/{sample}.clonotypes.TRA.txt" ),
        os.path.join( out_dir, "analyze/{sample}.clonotypes.TRB.txt" )
    run:
        sample_folder = os.path.join( out_dir, "analyze", wildcards.sample )
        # First build the command
        command = "mixcr analyze amplicon --species {} --starting-material {} --5-end {} --3-end {} --adapters {} --receptor-type {} {} {} {} {}".format(
            config["repertoire_assembly"]["species"],
            config["repertoire_assembly"]["starting_material"],
            config["repertoire_assembly"]["5-end"],
            config["repertoire_assembly"]["3-end"],
            config["repertoire_assembly"]["adapters"],
            config["repertoire_assembly"]["receptor_type"],
            config["repertoire_assembly"]["etc"],
            input[0],
            input[1],
            sample_folder
        )
        subprocess.call( command, shell=True )

rule remove_contamination:
    input:
        expand( os.path.join( out_dir, "umi_collapsed/{sample}_{read}.fastq.gz" ), out_dir=out_dir, sample=SAMPLES, read=["R1", "R2"] )
    params:
        input_folder = os.path.join( out_dir, "umi_collapsed" ),
        output_dir = os.path.join( out_dir, "cleaned" )
    output:
        expand( os.path.join( out_dir, "cleaned/{sample}_{read}.fastq.gz" ), out_dir=out_dir, sample=SAMPLES, read=["R1", "R2"] )
    shell:
        "python3 res/remove_contamination.py -t {config[remove_contamination][threshold]} {params.input_folder} {params.output_dir}"

rule umi_collapse_bulk:
    """ Assembles reads containing the same UMI.
    
    """

    input:
        expand( os.path.join( out_dir, "umi_extracted/{sample}_{read}.fastq.gz" ), out_dir=out_dir, sample=SAMPLES, read=["R1", "R2"] ),
        os.path.join( out_dir, "mig_statistics/estimates.txt" )
    params:
        umi_extraction_folder = os.path.join( out_dir, "umi_extracted/" ),
        mig_statistics_folder = os.path.join( out_dir, "mig_statistics/" ),
        output_folder = os.path.join( out_dir, "umi_collapsed/" )
    output:
        expand( os.path.join( out_dir, "umi_collapsed/{sample}_{read}.fastq.gz" ), out_dir=out_dir, sample=SAMPLES, read=["R1", "R2"] )
    run:
        command = "migec AssembleBatch -c {} {} {}".format( params.umi_extraction_folder, params.mig_statistics_folder, params.output_folder )
        subprocess.call( command, shell=True )

        #Remove the stupid options in the filename.
        command = 'for i in {}*.fastq.gz ; do mv "$i" "{}$(echo $i | rev | cut -f1 -d/ | rev | cut -f1,4,5 -d. )"; done'.format( params.output_folder, params.output_folder )
        subprocess.call( command, shell=True )

# Shelve this rule for now. MIGEC isn't functioning the way I'd like it too. A issue has been raised.
#rule umi_collapse:
#    input:
#        os.path.join( out_dir, "umi_extracted", "{sample}_R1.fastq.gz" ),
#        os.path.join( out_dir, "umi_extracted", "{sample}_R2.fastq.gz" ),
#        os.path.join( out_dir, "mig_statistics/estimates.txt" )
#    output:
#        os.path.join( out_dir, "umi_collapsed/{sample}_R1.fastq.gz" ),
#        os.path.join( out_dir, "umi_collapsed/{sample}_R2.fastq.gz" ),
#
#    params:
#        sample = "{sample}",
#        logfile = os.path.join( out_dir, "umi_collapsed/assemble.log.txt" )
#    run:
#        # Grab estimates from estimates.txt
#        estimates = pd.read_csv( input[2], sep="\t" )
#        overseq_threshold = int( estimates[estimates["#SAMPLE_ID"]==params.sample]["OVERSEQ_THRESHOLD"] )
#        quality_threshold = int( estimates[estimates["#SAMPLE_ID"]==params.sample]["UMI_QUAL_THRESHOLD"] )
#        # Call Assembly on build command
#        command = "migec Assemble -m {} -q {} --log-file {} --log-sample-name {} --log-sample-type paired {} -c {} {} {}".format( overseq_threshold,
#                                                                                                                                  quality_threshold,
#                                                                                                                                  params.logfile,
#                                                                                                                                  params.sample,
#                                                                                                                                  config["umi_collapse"]["assemble_options"],
#                                                                                                                                  input[0],
#                                                                                                                                  input[1],
#                                                                                                                                  os.path.join( out_dir, "umi_collapsed/" ) )
#        subprocess.call( command, shell=True )
#        # Oh and then because MIGEC attaches a stupid add tag to the file name, you've got to rename the files to make them work with snakemake.
#        command = "mv {} {}".format( os.path.join( out_dir, "umi_collapsed/", "{}_R1.t{}.fastq.gz".format( params.sample, overseq_threshold ) ), output[0] )
#        subprocess.call( command, shell=True )
#        command = "mv {} {}".format( os.path.join( out_dir, "umi_collapsed/", "{}_R2.t{}.fastq.gz".format( params.sample, overseq_threshold ) ), output[1] )
#        subprocess.call( command, shell=True )

rule mig_statistics:
    """ Generates molecular identifier group (MIG) size distribution statistics.
    
    Inputs:
    - UMI extracted fastqs  : folder containing all fastqs which have been processed with MIGEC Checkout
    - checkout.filelist.txt : file indicating the fastq locations of each sample. Used by MIGEC Histogram, etc.
    
    Outputs:
    - overseq.txt           : for each sample, lists the number of MIGs which have a given read coverage.
    - estimates.txt         : for each sample, lists the MIG size cutoff which dissects erroneous MIGs while retaining 
                            amplified ones
    - collisions.txt        : for each sample, lists the number of MIGs, which have a 1-mismatch UMI neighbour with a 
                            substantially higher count, which have a given read coverage.
    - pwm.txt               : position weight matrix (PWM) representation of all UMI sequences
    """

    input:
        expand( os.path.join( "{out_dir}", "umi_extracted", "{sample}_{read}.fastq.gz" ), out_dir=out_dir, sample=SAMPLES, read=["R1", "R2"] ),
        os.path.join( out_dir, "umi_extracted", "checkout.filelist.txt" )
    params:
        folder = os.path.join( out_dir, "umi_extracted" ),
        histogram_folder = os.path.join( out_dir, "mig_statistics")
    output:
        os.path.join( out_dir, "mig_statistics/estimates.txt" )
    shell:
        "migec Histogram {params.folder} {params.histogram_folder}"

rule cleanup_extraction:
    """ Collects umi extracted fastqs to a central location, generates files to trick MIGEC into thinking CheckoutBatch
    was performed, and removes superfluous files.
    
    Inputs:
    - Directory structure such that fastq pairs are found in folders with the name of their sample.
    
    Output:
    - Fastq pairs all in the same location
    - checkout.filelist.txt : file indicating the fastq locations of each sample. Used by MIGEC Histogram, etc.
    - checkout.log.txt      : Log file describing how many reads contained recognizable UMIs.
    
    Rationale:
    MIGEC doesn't allow you to specify a output prefix for any of its commands, so MIGEC Checkout will 
    override all of its generic outputs for each sample it processes. This is problematic because of the data loss, but 
    also because it disrupts parallel processing. To get around this, the pipeline creates a unique folder for each 
    sample processed, which need to be combined by this rule.
    """

    input:
        expand( os.path.join( "{out_dir}", "umi_extracted", "{sample}", "{sample}_{read}.fastq.gz" ), out_dir=out_dir, sample=SAMPLES, read=["R1", "R2"] )
    output:
        expand( os.path.join( "{out_dir}", "umi_extracted", "{sample}_{read}.fastq.gz" ), out_dir=out_dir, sample=SAMPLES, read=["R1", "R2"] ),
        os.path.join( out_dir, "umi_extracted", "checkout.filelist.txt" )
    run:
        # Move all fastq files from the individual sample folders to the parent folder
        command = "mv {} {}".format( os.path.join( out_dir, "umi_extracted", "*/*.fastq.gz" ), os.path.join( out_dir, "umi_extracted/" ) )
        subprocess.call( command, shell=True ) # THIS IS DANGEROUS, I guess...

        # Generate new filelist, everything else is superfluous.
        with open( os.path.join( out_dir, "umi_extracted", "checkout.filelist.txt" ), "w" ) as filelist:
            for i in SAMPLES:
                filelist_line = [i, "paired", os.path.join( out_dir, "umi_extracted", "{}_R1.fastq.gz".format( i ) ), os.path.join( out_dir, "umi_extracted", "{}_R2.fastq.gz".format( i ) ), "\n"]
                filelist.write( "\t".join( filelist_line ) )

        # Remove stupid undefined fastqs that MIGEC isn't supposed to keep.
        for i in ["undef-s_R2.fastq.gz", "undef-s_R1.fastq.gz", "undef-m_R2.fastq.gz", "undef-m_R1.fastq.gz"]:
            subprocess.call( ["rm", os.path.join( out_dir, "umi_extracted", i )] )
rule umi_extraction:
    """Iterates through the reads in a pair of fastq files and extracts the unique molecular identifier.
    
    Inputs:
    - fastq1 & fastq2   : the fastq files containing the reads for a single TCR library, i.e. after demultiplexing. 
                        Found automatically by snakemake
    - barcode           : A barcode sequences where upper and lower case letters mark seed and fuzzy-search region 
                        parts, respectively, and N characters mark UMI region to be extracted. Specified in the snakemake config file.
    - migec_options     : Optional parameters to pass through to the MIGEC Checkout command. Specified in snakemake 
                        config file. Defaults to -cute.
    - ouput_dir         : Folder to deposit output of MIGEC Checkout. 
    
    Outputs:
    - fastq1 & fastq2   : Depends migec_options, but by default output are input fastq files where each read where a UMI
    was successfully identified has its specific UMI added to the header. Unsuccessful reads are removed.
    
    Rationale: 
    There are two ways one can do this, in batch, or in manual steps. I've choosen to do it in manual step despite the 
    cons associated with it. For instance, manual processing forces creation of the barcode.txt file and then collection
    of the log file for each sample. The advantage, which makes it all worth it, is that by splitting processing up, 
    you can run each rule on its own job on Garibaldi. So this should be much faster. 
    """

    input:
        lambda wildcards: os.path.join( "{in_dir}".format( in_dir = in_dir ), SAMPLES[wildcards.sample][0] ),
        lambda wildcards: os.path.join( "{in_dir}".format( in_dir = in_dir ), SAMPLES[wildcards.sample][1] )
    output:
        os.path.join( "{out_dir}", "umi_extracted", "{sample}", "{sample}_R1.fastq.gz" ),
        os.path.join( "{out_dir}", "umi_extracted", "{sample}", "{sample}_R2.fastq.gz" )
    run:
        # Create the barcode.txt file which is just a tsv with the sample name, and barcode sequence. Needs to be placed
        # in the samples folder to avoid bullshit.
        barcode_loc = os.path.join( out_dir, "umi_extracted", wildcards.sample, "barcodes.txt" )
        with open( barcode_loc, "w" ) as barcodes_file:
            barcodes_string = [wildcards.sample, config["umi_extraction"]["barcode"], "", "\n"]
            barcodes_file.write( "\t".join( barcodes_string ) )

        # Call MIGEC to extract the UMIs using the barcode.txt file just created.
        subprocess.call( ["migec", "Checkout", config["umi_extraction"]["migec_options"], barcode_loc, input[0], input[1], os.path.join( out_dir, "umi_extracted", wildcards.sample )] )

