import json


disease = "17_11_2021_CDLeuven"
input_file = "pipeline_results_11_04_2024/CD/GO_terms_strange_all_CD_healthy_14_04_2024.tsv"
network_file = f"pipeline_results_11_04_2024/CD/{disease}_network_all_snps_13_03_2024.tsv"
master_table = f"files/iSNP_master_table_without_miRNA_{disease}_without_SNP_duplication.tsv"
go_file = "new_files/go_27_01_2022.json"
go_annotation_file = "new_files/go_annotations_27_01_2022.json"

healthy_specific_GOs = []
disease_specific_GOs = []
snps_dictionary = {}
goterms_snps_dictionary = {}
patients_dictionary = {}
final_dictionary = {}
patients_with_most_SNPs = []
go_dictionary = {}
go_annotations_dictionary = {}

with open(go_file, 'r') as go:

    for go_line in go:
        go_line = json.loads(go_line.strip())

        if go_line["name"] not in go_dictionary:
            go_dictionary[go_line["name"]] = go_line["id"]

with open(go_annotation_file, 'r') as go_annotation:

    for go_annotation_line in go_annotation:
        go_annotation_line = json.loads(go_annotation_line.strip())
        go_id = int(go_annotation_line["go_id"].split(":")[1])

        if go_id not in go_annotations_dictionary:
            go_annotations_dictionary[go_id] = []

        if go_annotation_line["id"] not in go_annotations_dictionary[go_id]:
            go_annotations_dictionary[go_id].append(go_annotation_line["id"])

with open(input_file, 'r') as i:

    i.readline()

    for line in i:
        line = line.strip().split('\t')

        if line[0] not in healthy_specific_GOs:
            line[0] = line[0].replace("_", " ")
            healthy_specific_GOs.append(line[0])

        if line[2] not in disease_specific_GOs:
            line[2] = line[2].replace("_", " ")
            disease_specific_GOs.append(line[2])

with open(network_file, 'r') as net:

    net.readline()

    for netline in net:
        netline = netline.strip().split('\t')
        source_go = netline[0].split("_")[0]
        target_go = netline[1].split("_")[0]
        snps = netline[3].split(",")

        if source_go in healthy_specific_GOs and target_go in disease_specific_GOs:
            interaction = f"{source_go}|||{target_go}"

            new_info = f"{source_go}|||{target_go}"

            if new_info not in goterms_snps_dictionary:
                goterms_snps_dictionary[new_info] = snps
            
            for snp in snps:
                snp = snp.upper()

                if snp not in snps_dictionary:
                    snps_dictionary[snp] = []

                snps_dictionary[snp].append(interaction)

with open(master_table, 'r') as master:

    for masterline in master:

        masterline = masterline.strip().split('\t')
        actual_snp = masterline[2]

        if actual_snp not in snps_dictionary:
            continue

        for x in range(5, len(masterline)):
            actual_patient_number = x - 4
            actual_patient = f"patient_{actual_patient_number}"

            if actual_patient not in patients_dictionary:
                patients_dictionary[actual_patient] = []

            if masterline[x] == "1":
                patients_dictionary[actual_patient].append(actual_snp)

# for k, v in patients_dictionary.items():

#     if len(v) == 6:
#         print(k, ": ", len(v), v)

#         if len(patients_with_most_SNPs) < 5:
#             patients_with_most_SNPs.append(k)

# print(patients_with_most_SNPs)

for pat, snplist in patients_dictionary.items():

    if pat not in final_dictionary:
        final_dictionary[pat] = []

    for actualsnpvariant in snplist:

        for actualgoterm in snps_dictionary[actualsnpvariant]:

            if actualgoterm not in final_dictionary[pat]:
                final_dictionary[pat].append(actualgoterm)

for k, v in final_dictionary.items():
    print(k)

    # if k in patients_with_most_SNPs:
    #     print(k)
    
    # # if k == "patient_40":
    # if k == "patient_35":

    with open(f"pipeline_results_11_04_2024/CD/cytoscape_networks/CD_{k}.tsv", 'w') as out:

        out.write("HealthyGOTerm\tDiseaseGOTerm\tHealthyGOTermTFs\tDiseaseGOTermTFs\tRelevantSNPs")

        for int in v:
            final_source = int.split("|||")[0]
            final_target = int.split("|||")[1]

            final_source_go_id = go_dictionary[final_source]
            final_target_go_id = go_dictionary[final_target]
            final_source_tfs = go_annotations_dictionary[final_source_go_id]
            final_target_tfs = go_annotations_dictionary[final_target_go_id]

            out.write(final_source + '\t' + final_target + '\t' + ",".join(final_source_tfs) + '\t' + ",".join(final_target_tfs) + '\t' + ",".join(goterms_snps_dictionary[int]) + '\n')
