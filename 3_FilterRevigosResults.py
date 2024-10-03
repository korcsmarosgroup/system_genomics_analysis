# A script to filter the Revigo results, according to the IBD-related GO term list (created by Dezso & Tamas)

import os

input_folder = "pipeline_results_11_04_2024/CD/Revigo_results/"
output_folder = "pipeline_results_11_04_2024/CD/Revigo_results_notfiltered/"
filtered_IBD_related_go_terms = "files/GO_terms_after_filtering_2022_05_10.txt"

# filtered_GO_terms = []

# with open(filtered_IBD_related_go_terms, 'r') as filtered:

#     for filtered_line in filtered:
#         go_term_name = filtered_line.strip()

#         if go_term_name not in filtered_GO_terms:
#             filtered_GO_terms.append(go_term_name)

for snp_folder in os.listdir(input_folder):

    if "." in snp_folder:
        continue

    output_file = os.path.join(output_folder, f"{snp_folder}.tsv")
    healthy_go_term_names = []
    non_healthy_go_term_names = []

    for filename in os.listdir(os.path.join(input_folder, snp_folder)):

        if filename != "healthy_processed.tsv" and filename != "non_healthy_processed.tsv":
            continue

        filepath = os.path.join(input_folder, snp_folder, filename)
        condition = "healthy"

        if filename == "non_healthy_processed.tsv":
            condition = "non_healthy"

        with open(filepath, 'r') as i:
            i.readline()

            for line in i:
                new_line = line.strip()
                i_line = new_line.split('\t')

                # if i_line[0] in filtered_GO_terms:

                if condition == "healthy":

                    if i_line[0] not in healthy_go_term_names:
                        healthy_go_term_names.append(i_line[0])

                if condition == "non_healthy":

                    if i_line[0] not in non_healthy_go_term_names:
                        non_healthy_go_term_names.append(i_line[0])

    if len(healthy_go_term_names) == 0 or len(non_healthy_go_term_names) == 0:
        continue

    with open(output_file, 'w') as out:

        out.write("Healthy.GO.Terms" + '\t' + "Non-healthy.GO.Terms" + '\n')
        go_terms_max_number = len(healthy_go_term_names)
        print(snp_folder, "Healthy: ", len(healthy_go_term_names), "Unhealthy: ", len(non_healthy_go_term_names))

        if len(non_healthy_go_term_names) > len(healthy_go_term_names):
            go_terms_max_number = len(non_healthy_go_term_names)

        if len(healthy_go_term_names) != 0 and len(non_healthy_go_term_names) != 0:

            for index in range(0, go_terms_max_number):

                if index < len(healthy_go_term_names) and index < len(non_healthy_go_term_names):
                    out.write(healthy_go_term_names[index] + '\t' + non_healthy_go_term_names[index] + '\n')

                if index < len(healthy_go_term_names) and index >= len(non_healthy_go_term_names):
                    out.write(healthy_go_term_names[index] + '\t' + "-" + '\n')

                if index >= len(healthy_go_term_names) and index < len(non_healthy_go_term_names):
                    out.write("-" + '\t' + non_healthy_go_term_names[index] + '\n')

        if len(healthy_go_term_names) != 0 and len(non_healthy_go_term_names) == 0:

            for index in range(0, len(healthy_go_term_names)):
                out.write(healthy_go_term_names[index] + '\t' + "-" + '\n')

        if len(healthy_go_term_names) == 0 and len(non_healthy_go_term_names) != 0:

            for index in range(0, len(non_healthy_go_term_names)):
                out.write("-" + '\t' + non_healthy_go_term_names[index] + '\n')
