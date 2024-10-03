# A script, which take the output of the previous script and put together the target genes per SNP
# Input file: the output of the *CollectTFSwitches.py* script
# Output file: a processed file, where the SNPs contains their target genes in one row

# Import packages
import argparse
import logging

# Create log level
logging.basicConfig(format = '%(asctime)s | %(levelname)s: %(message)s', level = logging.INFO)

# Create argparse argumentum for the terminal parameters
parser = argparse.ArgumentParser()
parser.add_argument('tf_switches_file')
parser.add_argument('processed_tf_switches_file')
argument = parser.parse_args()

# Create a dictionary, where collect all of the target genes per SNP
SNP_target_genes = {}

# Create a check array for the SNP to avoid duplicates
SNP_check_array = []

# Open the TF switches file
with open(argument.tf_switches_file, 'r') as i:

    # Skip the header
    i.readline()

    # Iterate through all of the line in the TF switches file
    for i_line in i:
        i_line = i_line.strip().split('\t')

        # Define the necessaray values
        snp_id = i_line[0]
        target_gene = i_line[1]

        # Check if the SNP is in the SNP target genes dictionary
        if snp_id not in SNP_target_genes:

            # Create an array for the target genes
            SNP_target_genes[snp_id] = []

        # Check if the target gene is in the array
        if target_gene not in SNP_target_genes[snp_id]:

            # Then put it into it
            SNP_target_genes[snp_id].append(target_gene)

# Open again the TF switches file and the output as well
with open(argument.tf_switches_file, 'r') as ii, open(argument.processed_tf_switches_file, 'w') as out:

    # Skip the header
    ii.readline()

    # Iterate through all of the line in the TF switches file
    for ii_line in ii:
        ii_line = ii_line.strip().split('\t')

        # Define the necessaray values
        snp_id = ii_line[0]
        target_gene = ii_line[1]

        # Check if the SNP is in the SNP check array
        if snp_id not in SNP_check_array:

            # Then put it into it
            SNP_check_array.append(snp_id)

            # Write out the line with all the target genes of the SNP
            out.write(snp_id + '\t' + ",".join(SNP_target_genes[snp_id]) + '\t' + ii_line[2] + '\t' + ii_line[3] + '\n')
