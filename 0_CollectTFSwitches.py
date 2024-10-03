# A script to collect all of the TF switches from the master tables from the Leuven cohorts
# Input file: A master table file without miRNAs from a cohort
# Output file: File, which will contain all of the TF switches between healthy and non-healthy conditions
# Uniprot IDs file: A file, which will contain all of the TF uniprot IDs from the master table
# Mapping file: A file, which contains the mapping information between uniprot and gene names

#disease = "05_11_2021_UCLeuven"
disease = "17_11_2021_CDLeuven"

# Define input, output and other files
input_file = f"files/iSNP_master_table_without_miRNA_{disease}.txt"
output_file = f"pipeline_results_11_04_2024/CD/TF_switches_results_{disease}.tsv"
uniprot_ids_file = f"files/uniprot_ids_for_mapping_{disease}.txt"
mapping_file = f"files/uniprot_to_genename_mapping_{disease}.tsv"

# Define dictionaries for wild type and for mutant type, where I store all of the SNP_targetgene - TF interactions
# Format: wild/mutant: {SNP_targetgene: [TF1, TF2, TF3]}
wild_dictionary = {}
mutant_dictionary = {}

# Define the difference dictionary
# Format: SNP_targetgene: {wild: [TF1, TF2, TF3], mutant: [TF4, TF5, TF6]}
differences = {}

# Define the mapping dictionary
# Format: Uniprot ID: {gene name}
mapping_dictionary = {}

# Array to collect all of the TF and target gene uniprot IDs
uniprot_ids = []

def FillUpDictionary(tfid, targetgeneid, snpid, dictionary, uniarray):
    """
    Function to fill up the wild type and the mutant type dictionaries
    """

    # Define the wild type SNP - target gene connection
    snp_targetgene = f"{snpid}_{targetgeneid}"

    # Check if the SNP - target gene connection is in the wild type dictionary
    if snp_targetgene not in dictionary:

        # If not, then put it there, and create an empty array for it to append the TF's later
        dictionary[snp_targetgene] = []

    # Check if the TF is inside the SNP - target gene array in the wild type dictionary
    if tfid not in dictionary[snp_targetgene]:

        # If not, then put the TF into the array
        dictionary[snp_targetgene].append(tfid)

    # Check if the TF is in the uniprot IDs array or not, if not, then put it into it
    if tfid not in uniarray:
        uniarray.append(tfid)

    # Check if the target gene is in the uniprot IDs array or not, if not, then put it into it
    if targetgeneid not in uniarray:
        uniarray.append(targetgeneid)


def SearchDifferences(difference_dictionary, first_dictionary, second_dictionary, condition):
    """
    Function to search all of the TFBS differences between the 2 conditions (healthy vs. mutant)
    """

    # Iterate through all of the SNP - target gene connections in the first dictionary
    for snptargetgene, tfs in first_dictionary.items():

        # If the SNP - target gene connection is not in the second dictionary
        if snptargetgene not in second_dictionary:

            # Check if the SNP - target gene connection is in the difference dictionary or not
            # If not, then create a key for it
            if snptargetgene not in difference_dictionary:
                difference_dictionary[snptargetgene] = {"wild": [], "mutant": []}

            # Iterate through the TF's inside the SNP - target gene connection
            for t in tfs:

                # Check if the TF is in the SNP - target gene connection of the difference dictionary
                # If not, then put it into it
                if t not in difference_dictionary[snptargetgene][condition]:
                    difference_dictionary[snptargetgene][condition].append(t)

        # If the SNP - target gene connection is in the second dictionary
        else:

            # Check if the SNP - target gene connection is in the difference dictionary or not
            # If not, then create a key for it
            if snptargetgene not in difference_dictionary:
                difference_dictionary[snptargetgene] = {"wild": [], "mutant": []}

            # Iterate through the TF's inside the SNP - target gene connection
            for tt in tfs:

                # Check if the TF is in the SNP - target gene connection of the second dictionary
                # If not, then put it into it
                if tt not in second_dictionary[snptargetgene]:
                    difference_dictionary[snptargetgene][condition].append(tt)


def CollectDifferences(difference_dictionary, condition, condition_array):
    """
    Function to collect all of the TFBS differences between the 2 conditions (healthy vs. mutant) per SNP per gene
    """

    # If the array is not empty
    if len(difference_dictionary[snptargetgene][condition]) != 0:

        # Iterate through the TFs in the wild condition, then put the TF name (using the mapping dictionary)
        # into the wild_array
        for condition_tf in difference_dictionary[snptargetgene][condition]:
                condition_array.append(mapping_dictionary[condition_tf])

    # If the array is empty, then put a single "-" into it
    if len(difference_dictionary[snptargetgene][condition]) == 0:
        condition_array.append("-")


# Open the input file and the file, where I collect all of the TF and target gene uniprot IDs for mapping
with open(input_file, 'r') as i, open(uniprot_ids_file, 'w') as tfuni:

    # Iterate through all of the lines in the input file
    for line in i:
        line = line.strip().split('\t')

        # Define the necessary values (TF ID, target gene ID, SNP ID, condition - healhty or mutant -)
        tf = line[0]
        target_gene = line[1]
        snp_id = line[2]
        condition = line[3]

        # Check the first column uniprot ID
        if "/" in tf:
            continue 

        # Fill up the wild type dictionary
        if condition == "WT":
            FillUpDictionary(tf, target_gene, snp_id, wild_dictionary, uniprot_ids)

        # Fill up the mutant type dictionary
        if condition == "MUT":
            FillUpDictionary(tf, target_gene, snp_id, mutant_dictionary, uniprot_ids)

    # Write out all of the uniprot IDs for mapping
    for uni_id in uniprot_ids:
        tfuni.write(uni_id + '\n')

# Search the difference from the healthy condition
SearchDifferences(differences, wild_dictionary, mutant_dictionary, "wild")

# Search the difference from the mutant condition
SearchDifferences(differences, mutant_dictionary, wild_dictionary, "mutant")

# Make a dictionary from the mapping file, open the mapping file
with open(mapping_file, 'r') as mapping:

    # Skip the header
    mapping.readline()

    # Iterate through all of the lines in the mapping file
    for mapping_line in mapping:
        mapping_line = mapping_line.strip().split('\t')

        # Define the necessary values (uniprot id, genename)
        u_id = mapping_line[0]
        genename = mapping_line[1]

        # Check if the uniprot id is in the mapping dictionary or not, if not then put it in with the genename
        if u_id not in mapping_dictionary:
            mapping_dictionary[u_id] = genename

# Write out the results into the output file
with open(output_file, "w") as out:

    # Write out the header
    out.write("SNP.ID" + '\t' + "Target.Gene.Name" + '\t' + "Healthy.TFs" + '\t' + "UC.TFs" + '\n')

    # Iterate through the differences dictionary
    for snptargetgene in differences:

        # Define the necessary values (SNP ID, target gene uniprot ID)
        snpid = snptargetgene.split("_")[0]
        targetgeneid = snptargetgene.split("_")[1]

        # Define the arrays, where I can collect the wild and mutant condition TFs
        wild_array = []
        mutant_array = []

        # Collect all of the wild condition TFs
        CollectDifferences(differences, "wild", wild_array)

        # Collect all of the mutant condition TFs
        CollectDifferences(differences, "mutant", mutant_array)

        # Write out the difference into the output file
        out.write(snpid + '\t' + mapping_dictionary[targetgeneid] + '\t' + ",".join(wild_array) + '\t' + 
                    ",".join(mutant_array) + '\n')
