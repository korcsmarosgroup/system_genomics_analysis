import sys


#disease = "05_11_2021_UCLeuven"
disease = "17_11_2021_CDLeuven"
input_file = "pipeline_results_11_04_2024/CD/CD_GO_switches_filtered_tamma_13_04_2024.tsv"
master_file = f"files/iSNP_master_table_without_miRNA_{disease}_without_SNP_duplication.tsv"
output_file = f"pipeline_results_11_04_2024/CD/{disease}_network_all_SNP_associated_clustering_input_13_04_2024.tsv"
heatmap_input = f"pipeline_results_11_04_2024/CD/{disease}_network_all_SNP_associated_heatmap_input_13_04_2024.tsv"
lines_in_input = sum(1 for line in open(input_file))
patient_dictionary = {}

with open(master_file, 'r') as master:

    for master_line in master:
        master_line = master_line.strip().split('\t')
        actual_snp = master_line[2].lower()

        if actual_snp not in patient_dictionary:
            patient_dictionary[actual_snp] = []

        for x in range(5, len(master_line)):
            patient_dictionary[actual_snp].append(master_line[x])

# for k, v in patient_dictionary.items():

#     print(k, " : ", len(v))

index = 1
with open(input_file, 'r') as i, open(output_file, 'w') as out, open(heatmap_input, 'w') as heatmap:

    i.readline()
    patients_ids = []
    for z in range(1, 1696): # For CD
    #for z in range(1, 942): # For UC
        patients_ids.append(f"Patient_{z}")
    out.write("Source|Target" + '\t' + '\t'.join(patients_ids) + '\n')
    heatmap.write("GoTerms" + '\t' + "EffectingSNP" + '\n')

    for line in i:
        line = line.strip().split('\t')
        source = line[0]
        target = line[1]
        snps = line[3].split(",")
        patient_info = {}

        new_source = source.replace(" ", "_")
        new_target = target.replace(" ", "_")

        for snp in snps:
            patient_info[snp] = patient_dictionary[snp]

        output_patient_info = []

        if len(snps) == 1:
            heatmap_value = snps[0]

            for y in range(0, len(patient_info[snps[0]])):
                first_number = patient_info[snps[0]][y]
                output_patient_info.append(first_number)

        elif len(snps) == 2:
            heatmap_value = "Multiple"

            for y in range(0, len(patient_info[snps[0]])):
                first_number = patient_info[snps[0]][y]
                second_number = patient_info[snps[1]][y]

                if first_number == second_number:
                    output_patient_info.append(first_number)
                elif first_number != second_number:
                    output_patient_info.append(str(1))

        elif len(snps) == 3:
            heatmap_value = "Multiple"

            for y in range(0, len(patient_info[snps[0]])):
                first_number = patient_info[snps[0]][y]
                second_number = patient_info[snps[1]][y]
                third_number = patient_info[snps[2]][y]

                if first_number == second_number and second_number == third_number:
                    output_patient_info.append(first_number)
                elif first_number != second_number or first_number != third_number:
                    output_patient_info.append(str(1))

        out.write(f"{new_source}|{new_target}" + '\t' + '\t'.join(output_patient_info) + '\n')
        heatmap.write(f"{new_source}|{new_target}" + '\t' + heatmap_value + '\n')

        percentage = int(index / lines_in_input * 100)
        sys.stdout.write(f"\r{index} / {lines_in_input} ---> {percentage}%")
        index += 1

sys.stdout.write(f"\r{lines_in_input} / {lines_in_input} ---> 100%")
print("\nDone!")
