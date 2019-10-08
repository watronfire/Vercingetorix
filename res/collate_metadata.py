import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument( "input_dir", help="Pipeline output directory" )
parser.add_argument( "output_prefix", help="Name for output csv file" )

args = parser.parse_args()
input_dir = args.input_dir
output_file = args.output_prefix

samples = [_ for _ in os.listdir( os.path.join( input_dir, "umi_extracted" ) ) if os.path.isdir( os.path.join( input_dir, "umi_extracted", _ ) )]

extraction_list = list()

for _ in samples:
    temp = pd.read_csv( os.path.join( input_dir, "umi_extracted", _, "checkout.log.txt" ), sep="\t" )
    umis_identified = int( temp.loc[temp["SAMPLE"]==_, "MASTER"] )
    total_reads = umis_identified + int( temp.loc[temp["SAMPLE"]=="undef-m", "MASTER"] )
    extraction_list.append( [_, total_reads, umis_identified ] )
metadata_df = pd.DataFrame( extraction_list, columns=["Sample", "Total_Reads", "Barcode_Reads" ] )
metadata_df.set_index( "Sample", inplace=True )

ms = pd.read_csv( os.path.join( input_dir, "mig_statistics", "estimates.txt"), sep="\t" )
ms.rename( columns={"#SAMPLE_ID" : "Sample",
                    "SAMPLE_TYPE" : "Sample_Type",
                    "TOTAL_MIGS" : "Total_MIGS",
                    "OVERSEQ_THRESHOLD" : "Overseq_Threshold",
                    "COLLISION_THRESHOLD" : "Collision_Threshold",
                    "UMI_QUAL_THRESHOLD" : "UMI_Quality_Threshold",
                    "UMI_LEN" : "UMI_Length" }, inplace=True )
ms.drop( "TOTAL_READS", inplace=True, axis=1 )
ms.set_index( "Sample", inplace=True )
metadata_df = pd.concat( [metadata_df, ms], axis=1 )

uc = pd.read_csv( os.path.join( input_dir, "umi_collapsed", "assemble.log.txt" ), sep="\t" )
uc.drop( ["SAMPLE_TYPE", "INPUT_FASTQ1", 'MIGS_DROPPED_OVERSEQ_2', 'READS_DROPPED_OVERSEQ_2', 'READS_DROPPED_COLLISION_2', 'MIGS_DROPPED_COLLISION_2', "INPUT_FASTQ2", "OUTPUT_ASSEMBLY1", "OUTPUT_ASSEMBLY2", "MIG_COUNT_THRESHOLD", "MIGS_TOTAL"], axis=1, inplace=True )
uc.rename( columns={ '#SAMPLE_ID' : 'Sample',
                     'MIGS_GOOD_FASTQ1' : 'MIGs_Good_R1',
                     'MIGS_GOOD_FASTQ2' : 'MIGs_Good_R2',
                     'MIGS_GOOD_TOTAL' : 'MIGs_Good_Total',
                     'READS_GOOD_FASTQ1' : 'Reads_Good_R1',
                     'READS_GOOD_FASTQ2' : 'Reads_Good_R2',
                     'READS_GOOD_TOTAL' : 'Reads_Good_Total',
                     'READS_TOTAL' : 'Total_Reads_Collapsed',
                     'READS_DROPPED_WITHIN_MIG_1' : 'Reads_Dropped_Within_MIG_1',
                     'READS_DROPPED_WITHIN_MIG_2' : 'Reads_Dropped_Within_MIG_2',
                     'MIGS_DROPPED_OVERSEQ_1' : 'MIGs_Dropped_Overseq',
                     'READS_DROPPED_OVERSEQ_1' : 'Reads_Dropped_Overseq',
                     'MIGS_DROPPED_COLLISION_1' : 'MIGs_Dropped_Collision',
                     'READS_DROPPED_COLLISION_1' : 'Reads_Dropped_Collision' }, inplace=True )
uc.set_index( "Sample", inplace=True )
metadata_df = pd.concat( [metadata_df, uc], axis=1 )

cr = pd.read_csv( os.path.join( input_dir, "cleaned/contamination.log.txt" ) )
cr.drop( ["TOTAL_READS"], axis=1, inplace=True )
cr.rename( columns={ "SAMPLE_ID" : "Sample",
                   "CONTAMINATED_READS" : "Contaminated_MIGs",
                   "GOOD_READS" : "Good_MIGs",
                   "GOOD_READS_PCT" : "Good_Reads_PCT" }, inplace=True )
cr.set_index( "Sample", inplace=True )
metadata_df = pd.concat( [metadata_df, cr ], axis=1, sort=True )

df = list()
for sample in samples:
    with open( os.path.join( input_dir, "analyze", "{}.report".format( sample ) ), "r" ) as ar:
        value_dict = dict()
        for _ in ar:
            line_split = _.split(":")
            if len( line_split ) > 1:
                factor = line_split[0].split( "," )[0].replace( " ", "_" )
                value_split = line_split[1].split( "(" )
                value =  value_split[0].strip()
                try:
                    if "." in value:
                        value_dict[factor] = float( value )
                    else:
                        value_dict[factor] = int( value )
                except ValueError:
                    continue
        df.append( pd.DataFrame( value_dict, index=[sample] ) )

df = pd.concat( df, sort=True )
df.drop( ["Total_sequencing_reads"], axis=1, inplace=True )
metadata_df = pd.concat( [metadata_df, df], axis=1, sort=True)

metadata_df.to_csv( output_file, sep="," )