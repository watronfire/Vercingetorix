import gzip
from Bio import SeqIO
import os
from fnmatch import fnmatch
from collections import OrderedDict
import Levenshtein
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument( "input_dir", help="directory where contaminated fastq file pairs are located" )
parser.add_argument( "output_dir", help="directory where cleaned fastq file pairs will be deposited" )
parser.add_argument( "-l", "--log_prefix", help="prefix to attach to logfiles" )
parser.add_argument( "-t", "--threshold", help="maximum distance allowed within a molecular identifier group (Default: 40)", type=int )

args = parser.parse_args()

in_dir = args.input_dir
out_dir = args.output_dir
cmdfile_loc = args.log_prefix + ".cmd.txt" if args.log_prefix else os.path.join( args.output_dir, "contamination.cmd.txt" )
logfile_loc = args.log_prefix + ".log.txt" if args.log_prefix else os.path.join( args.output_dir, "contamination.log.txt" )
threshold = args.threshold if args.threshold else 40

# Write the command to file.
with open( os.path.join( out_dir, cmdfile_loc ), "w" ) as cmdfile:
    cmdfile.write( " ".join( sys.argv ) )

SAMPLES = OrderedDict()

# Identified fastq pairs
for f in os.listdir( in_dir ):
    if fnmatch( f, "*.fastq.gz" ):
        fsplit = f.split( "_" )
        if fsplit[0] not in SAMPLES:
            SAMPLES[fsplit[0]] = [f]
        else:
            SAMPLES[fsplit[0]].append( f )
# Order file list for R1 is in the zeroeth index
for _ in SAMPLES:
    SAMPLES[_].sort()

observed_umis = OrderedDict()
contaminant_reads = list()
# Iterate through all the samples.
print( "Building UMI database for {} samples...".format( len( SAMPLES ) ) )
for _ in SAMPLES.items():
    with gzip.open( os.path.join( in_dir, _[1][1] ), "rt" ) as handle:
        # Iterate through all the reads in the sample
        for record in SeqIO.parse(handle, "fastq"):
            # Identify the UMI
            umi = record.description.split( ":" )[1]
            # check if this UMI has been observed before.
            if umi in observed_umis:
                # If it has then check the genetic distance between sequences associated with that UMI.
                for j in observed_umis[umi]:
                    if Levenshtein.distance( str( record.seq ), j[0] ) < 40:
                        # If the genetic distance is less than 40, likely result from the same molecule, i.e. contamination.
                        # check if UMI is in contaminant reads yet.
                        contaminant_reads.append( "{}-{}".format( _[0], record.id ) )
                        contaminant_reads.append( j[1] )
                    else:
                        observed_umis[umi].append( ( str( record.seq ), "{}-{}".format( _[0], record.id ) ) )
            else:
                observed_umis[umi] = [( str( record.seq ), "{}-{}".format( _[0], record.id ) )]

print( "UMI database completed." )
print( "{} UMIs observed.".format( len( observed_umis ) ) )
contaminations = len( list( set( contaminant_reads ) ) )
# This looks a little high but I'm not sure how else to check it.
print( "{} ({:.2f}%) contaminated sequences found.".format( contaminations, 100 * contaminations / len( observed_umis ) ) )

log_string = list()

for _ in SAMPLES.items():
    with gzip.open( os.path.join( in_dir, SAMPLES[_[0]][0] ), "rt" ) as handleR1, \
            gzip.open( os.path.join( in_dir, SAMPLES[_[0]][1] ), "rt" ) as handleR2, \
            gzip.open( os.path.join( out_dir, SAMPLES[_[0]][0] ), "wt" ) as outputR1, \
            gzip.open( os.path.join( out_dir, SAMPLES[_[0]][1] ), "wt" ) as outputR2:
        count = 0
        removed_reads = 0
        for records in zip( SeqIO.parse( handleR1, "fastq" ), SeqIO.parse( handleR2, "fastq" ) ):
            if "{}-{}".format( _[0], records[0].id ) not in contaminant_reads:
                outputR1.write( records[0].format( "fastq" ) )
                outputR2.write( records[1].format( "fastq" ) )
            else:
                removed_reads += 1
            count += 1
    log_string.append( [_[0], str( count ), str( removed_reads ), str( count - removed_reads), "{:.2f}%".format( 100 * ( count - removed_reads ) / count )] )

with open( os.path.join( out_dir, logfile_loc ), "w" ) as logfile:
    logfile.write( ",".join([ "SAMPLE_ID", "TOTAL_READS", "CONTAMINATED_READS", "GOOD_READS", "GOOD_READS_PCT"] ) + "\n" )
    for i in log_string:
        logfile.write( ",".join( i ) + "\n" )