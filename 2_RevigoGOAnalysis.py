# A script to check GO differences between healthy and non-healthy conditions
# Master table file: The master table, which was created by Dezső's script and contains all of the TF-target gene interactions from the cohort of interest
# TF switches file: The output file from the *CollectTFSwitches.py* script
# Disease all TFs file: A file, which contains all of the TF's from the master table of interest -> background for the enrichment analysis
# Mapping file: Contains the mapping information between gene names and unirpot ID's

# Import packages
import subprocess
import argparse
import logging
import json
import sys
import os

# Create log level
logging.basicConfig(format = '%(asctime)s | %(levelname)s: %(message)s', level = logging.INFO)

# Create argparse argumentum for the terminal parameters
parser = argparse.ArgumentParser()
parser.add_argument('tf_switches_file')
parser.add_argument('mapping_file')
parser.add_argument('go_annotations_file')
argument = parser.parse_args()

# Create a dictionary for the gene name - uniprot ID mapping
mapping_dictionary = {}

# Create a dictionary for the GO annotations of the TFs
go_annotations = {}

# Create a checking array for the SNPs to avoid duplicate SNP row analysis
check_snps = []

# Define the number of the TF switches file's lines
lines_in_input = sum(1 for line in open(argument.tf_switches_file))


def CollectGOIDs(condition_specific_tf_list, mapping_dictionary, condition_specific_helper_array):
    """
    Function to collect all of the condition specific TFs GO IDs and put them for the helper array
    """

    # Iterate through the TF list
    for tf in condition_specific_tf_list:

        # Check if the TF is inside the mapping dictionary
        if tf in mapping_dictionary:

            # Define the TF uniprot ID
            TF_uniprot_ID = mapping_dictionary[tf]

            # Check if the TF uniprot ID inside the GO annotations dictionary
            if TF_uniprot_ID in go_annotations:

                # Iterate through all of the GO IDs of the TF
                for go_id in go_annotations[TF_uniprot_ID]:

                    # Check if the GO ID is in the condition specific helper array
                    if go_id not in condition_specific_helper_array:

                        # Then put it into it
                        condition_specific_helper_array.append(go_id)


def CalculateDifferences(first_helper_array, second_helper_array, condition_specific_difference_array):
    """
    Function to calculate the differences between the two condition specific helper arrays with GO IDs
    """

    # Iterate through the condition specific GO IDs
    for specific_go_id in first_helper_array:

        # Check if the GO ID is in the other condition specific helper array
        if specific_go_id not in second_helper_array:

            # Check if the GO ID is inside the differnce helper array
            if specific_go_id not in condition_specific_difference_array:

                # Check if the specific GO ID is in the condition specific difference array
                if specific_go_id not in condition_specific_difference_array:

                    # Then put it into the differences condition specific helper array
                    condition_specific_difference_array.append(specific_go_id)


def GetTargetGeneGOIDs(gene_go_helper_file, gene_names_list, mapping_dict, go_dict):
    """
    Function to collect all the GO IDs of the given target gene
    """

    # Open the helper file
    with open(gene_go_helper_file, 'w') as gene_go_helper:

        # Write out the header to the helper file
        gene_go_helper.write("GOID.TFs" + '\n')

        # Iterate through the gene names
        for genename in gene_names_list:

            # Get the Uniprot ID of the gene name
            genename_uniprot_id = mapping_dict[genename]

            # Check if the gene name uniprot ID inside the GO annotations dictionary
            if genename_uniprot_id in go_dict:

                # Iterate through the GO IDs of the gene name uniprot ID
                for genename_go in go_dict[genename_uniprot_id]:

                    # Write them out to the helper file
                    gene_go_helper.write(genename_go + '\n')


def WriteOutGOIDsTOHelperFile(helper_file, condition_specific_difference_array):
    """
    Function to write out all of the TF and condition specific GO IDs to the helper file 
    """

    # Open the helper file to write out the GO IDs of the TFs for the Revigo analysis
    with open(helper_file, 'w') as helper_file:

        # Write out the header to the helper file
        helper_file.write("GOID.TFs" + '\n')

        # Iterate through the condition specific GO IDs
        for difference_go_id in condition_specific_difference_array:

            # Write out the GO ID to the helper file
            helper_file.write(difference_go_id + '\n')


def ProcessTSVFile(condition_specific_tsv_file, parent_dictionary, children_dictionary):
    """"
    Function to process the TSV files after the Revigo analysis
    """

    # Open the TSV file after the Revigo analysis
    with open(condition_specific_tsv_file, 'r') as tsv_file:

        # Skip the header in the TSV file
        tsv_file.readline()

        # Iterate through the lines in the TSV file
        for tsv_file_line in tsv_file:
            tsv_file_line = tsv_file_line.strip()

            # Delete the '"' symbol from the lines
            new_tsv_file_line = tsv_file_line.replace('"', "")
            new_tsv_file_line = new_tsv_file_line.split('\t')

            # Define the necessary values
            cluster_id = new_tsv_file_line[2]
            cluster_parent = new_tsv_file_line[7]
            actual_go_term_name = new_tsv_file_line[6]

            # Check if the cluster ID is in the parent dictionary
            if cluster_id not in parent_dictionary:

                # Put the parent GO term name into cluster ID of the parent dictionary
                parent_dictionary[cluster_id] = cluster_parent

            # Check if the cluster ID is in the children dictionary
            if cluster_id not in children_dictionary:

                # Create an array for the children GO term names
                children_dictionary[cluster_id] = []

            # Check if the actual GO term name is in the cluster ID of the children dictionary
            if actual_go_term_name not in children_dictionary[cluster_id]:

                # Check if the actual GO term is a parent term
                if actual_go_term_name != cluster_parent:

                    # Then put it into the cluster ID of the children dictionary
                    children_dictionary[cluster_id].append(actual_go_term_name)


def WriteOutTheProcessedGOTerms(condition_specific_tsv_file, children_dictionary, condition_specific_output_file):
    """
    Function to write out the processed GO terms to a new output file
    """

    # Open the TSV file after the Revigo analysis
    with open(condition_specific_tsv_file, 'r') as tsv_file, open(condition_specific_output_file, 'w') as out:

        # Skip the header in the TSV file
        tsv_file.readline()

        # Write out the header to the output file
        out.write("GO.Term.Name" + '\t' + "GO.Term.ID" + '\t' + "GO.cluster" + '\t' + "GO.score" + '\t' + "Children.of.GO.Term" + '\n')

        # Iterate through the lines in the TSV file
        for tsv_file_line in tsv_file:
            tsv_file_line = tsv_file_line.strip()

            # Delete the '"' symbol from the lines
            new_tsv_file_line = tsv_file_line.replace('"', "")
            new_tsv_file_line = new_tsv_file_line.split('\t')

            # Define the necessary values
            cluster_id = new_tsv_file_line[2]
            cluster_parent = new_tsv_file_line[7]
            actual_go_term_name = new_tsv_file_line[6]

            # Check if the cluster ID is in the parent dictionary
            if actual_go_term_name == cluster_parent:

                # Write out the new line into the output + the children GO term names
                out.write(new_tsv_file_line[6] + '\t' + new_tsv_file_line[0] + '\t' + new_tsv_file_line[2] + '\t' + new_tsv_file_line[4] + '\t' + " | ".join(children_dictionary[cluster_id]) + '\n')


print("")
logging.info("Everything is fine, starting...")
logging.info("Make the mapping dictionary for genename to uniprot ID mapping")

# Make a dictionary from the mapping file, open the mapping file
# with open(argument.mapping_file, 'r') as mapping:

#     # Skip the header
#     mapping.readline()

#     # Iterate through all of the lines in the mapping file
#     for mapping_line in mapping:
#         mapping_line = json.loads(mapping_line.strip())

#         # Define the necessary values (uniprot id, genename)
#         if mapping_line["from_id_type"] == "genename_primary" and mapping_line["to_id_type"] == "uniprotac":
#             gene_name = mapping_line["from_id"].upper()
#             uniprotac = mapping_line["to_id"].upper()

#             # Check if the uniprot id is in the mapping dictionary or not, if not then put it in with the genename
#             if gene_name not in mapping_dictionary:
#                 mapping_dictionary[gene_name] = uniprotac

# logging.info("Create the GO annotations dictionary for the TFs")

# # Open the GO annotation file
# with open(argument.go_annotations_file, 'r') as go_annot:

#     # Iterate through all of the line in the GO annotation file
#     for go_annot_line in go_annot:
#         go_annot_line = json.loads(go_annot_line.strip())

#         # Check if the TF uniprot ID is in the GO annotation dictionary
#         if go_annot_line["id"] not in go_annotations:

#             # If not then create an array for it inside the dictionary
#             go_annotations[go_annot_line["id"]] = []

#         # Check if the GO ID is in the TF array of the GO annotation dictionary
#         if go_annot_line["go_id"] not in go_annotations[go_annot_line["id"]]:

#             # Then put it into the array
#             go_annotations[go_annot_line["id"]].append(go_annot_line["go_id"])

logging.info("Start the enrichment analysis and the Revigo analysis per SNP")

# Open the TF switches file and do the enrichment and revigo analysis per line
with open(argument.tf_switches_file, 'r') as tf_switches:

    # Skip the header
    # tf_switches.readline()

    # Define the values for follow the process
    index = 1

    # Iterate through all of the lines from the TF switches file
    for tf_switches_line in tf_switches:
        tf_switches_line = tf_switches_line.strip().split('\t')

        # Add one to the follower value
        percentage = int(index / lines_in_input * 100)
        sys.stdout.write(f"\r          {index} / {lines_in_input} lines are processed - {percentage}%")
        index += 1

        # Define the necessary values
        snp_id = tf_switches_line[0]

        # Create a folder of the SNP in the disease-specific folder
        snp_folder = f"pipeline_results_11_04_2024/CD/Revigo_results/{snp_id}"
        if not os.path.isdir(snp_folder):
            os.mkdir(snp_folder)

        target_gene_names_list = tf_switches_line[1].split(",")
        healthy_tfs_list = tf_switches_line[2].split(",")
        non_healhty_tfs_list = tf_switches_line[3].split(",")
        healthy_go_ids = []
        non_healthy_go_ids = []
        difference_healthy_go_ids = []
        difference_non_healthy_go_ids = []
        target_gene_helper_file = f"{snp_folder}/target_gene_go_ids.tsv"
        healthy_helper_file = f"{snp_folder}/healthy_go_ids.tsv"
        non_healthy_helper_file = f"{snp_folder}/non_healthy_go_ids.tsv"
        go_cluster_parent_healthy = {}
        go_cluster_children_healthy = {}
        go_cluster_parent_non_healthy = {}
        go_cluster_children_non_healthy = {}
        go_cluster_parent_target_genes = {}
        go_cluster_children_target_genes = {}

        # Check if the SNP was already processed or not
        if snp_id in check_snps:
            continue

        # Add SNP id to the checking array
        check_snps.append(snp_id)

        # Collect all of the target gene specific GO terms
        # GetTargetGeneGOIDs(target_gene_helper_file, target_gene_names_list, mapping_dictionary, go_annotations)

        # # Collect all of the healthy specific TFs GO IDs into helper array
        # CollectGOIDs(healthy_tfs_list, mapping_dictionary, healthy_go_ids)

        # # Collect all of the non-healthy specific TFs GO IDs into helper array
        # CollectGOIDs(non_healhty_tfs_list, mapping_dictionary, non_healthy_go_ids)

        # # Calculate the differences between the healthy and non-healthy helper arrays
        # CalculateDifferences(healthy_go_ids, non_healthy_go_ids, difference_healthy_go_ids)

        # # Calculate the differences between the non-healthy and healthy helper arrays
        # CalculateDifferences(non_healthy_go_ids, healthy_go_ids, difference_non_healthy_go_ids)

        # # Write out the condition healthy GO IDs to the helper file
        # WriteOutGOIDsTOHelperFile(healthy_helper_file, difference_healthy_go_ids)

        # # Write out the condition non-healthy GO IDs to the helper file
        # WriteOutGOIDsTOHelperFile(non_healthy_helper_file, difference_non_healthy_go_ids)

        # # Run the R script to translate IDs to Entrez and do the Revigo analysis
        # command = ["Rscript", "RevigoAnalysis.R", healthy_helper_file, non_healthy_helper_file, target_gene_helper_file, snp_folder]
        # p = subprocess.Popen(command, stderr = subprocess.PIPE, stdout = subprocess.PIPE)
        # my_stdout, my_stderr = p.communicate()

        # Process the SNP-specific healthy TSV file
        healthy_tsv_file_path = os.path.join(snp_folder, "healthy.tsv")

        # Check if the file exists
        if os.path.isfile(healthy_tsv_file_path):
            ProcessTSVFile(healthy_tsv_file_path, go_cluster_parent_healthy, go_cluster_children_healthy)

            # Write out the healthy processed GO terms into a new output file
            healthy_output_file_path = os.path.join(snp_folder, "healthy_processed.tsv")
            WriteOutTheProcessedGOTerms(healthy_tsv_file_path, go_cluster_children_healthy, healthy_output_file_path)

        # Process the SNP-specific non-healthy TSV file
        non_healthy_tsv_file_path = os.path.join(snp_folder, "non_healthy.tsv")

        # Check if the file exists
        if os.path.isfile(non_healthy_tsv_file_path):
            ProcessTSVFile(non_healthy_tsv_file_path, go_cluster_parent_non_healthy, go_cluster_children_non_healthy)

            # Write out the non-healthy processed GO terms into a new output file
            non_healthy_output_file_path = os.path.join(snp_folder, "non_healthy_processed.tsv")
            WriteOutTheProcessedGOTerms(non_healthy_tsv_file_path, go_cluster_children_non_healthy, non_healthy_output_file_path)

        # Process the SNP-specific target gene TSV file
        target_gene_tsv_file_path = os.path.join(snp_folder, "target_gene.tsv")

        # Check if the file exists
        if os.path.isfile(target_gene_tsv_file_path):
            ProcessTSVFile(target_gene_tsv_file_path, go_cluster_parent_target_genes, go_cluster_children_target_genes)

            # Write out the target genes processed GO terms into a new output file
            target_gene_output_file_path = os.path.join(snp_folder, "target_gene_processed.tsv")
            WriteOutTheProcessedGOTerms(target_gene_tsv_file_path, go_cluster_children_target_genes, target_gene_output_file_path)

    sys.stdout.write(f"\r          {lines_in_input} / {lines_in_input} lines are processed - 100%")
    print("")

logging.info("Finished, great job!")
print("")
