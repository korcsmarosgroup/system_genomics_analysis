# A script to create a network in cytoscape about the GO term connections
# Node color: 2 (healthy, non-healhty)
# Node size: How many TF has this GO term, in a disease-specific way
# Directed edges!
# Edge color: How many times is it appearing in the SNPs

import json
import os


#disease = "05_11_2021_UCLeuven"
disease = "17_11_2021_CDLeuven"
input_folder = f"pipeline_results_11_04_2024/CD/Revigo_results_notfiltered"
TF_switches_file = f"pipeline_results_11_04_2024/CD/TF_switches_results_{disease}_processed.tsv"
go_file = "new_files/go_27_01_2022.json"
go_annotations_file = "new_files/go_annotations_27_01_2022.json"
mapping_file = f"files/human_uniprot_mapping_version_2022_01.json"
output_file = f"pipeline_results_11_04_2024/CD/{disease}_network_all_snps_23_07_2024.tsv"
output_node_file = f"pipeline_results_11_04_2024/CD/{disease}_nodes_all_snps_23_07_2024.tsv"
mapping_dictionary = {}
go_dictionary = {}
go_annotations_dictionary = {}
tfs_go_terms = {}
interactions_dictionary = {}
all_go_terms = []
only_healthy_go_terms = []
only_nonhealthy_go_terms = []
go_term_numbers_in_TFs = {}
interactions_per_snp = {}

print("Make mapping dictionary")
with open(mapping_file, 'r') as mapping:

    mapping.readline()

    for mapping_line in mapping:
        mapping_line = json.loads(mapping_line.strip())

        if mapping_line["from_id_type"] == "genename_primary" and mapping_line["to_id_type"] == "uniprotac":
            gene_name = mapping_line["from_id"].upper()
            uniprotac = mapping_line["to_id"].upper()

            if gene_name not in mapping_dictionary:
                mapping_dictionary[gene_name] = uniprotac

# with open(mapping_file, 'r') as mapping:

#     mapping.readline()

#     for mapping_line in mapping:
#         mapping_line = mapping_line.strip().split('\t')
#         uni_id = mapping_line[0]
#         g_name = mapping_line[1]

#         if g_name not in mapping_dictionary:
#             mapping_dictionary[g_name] = uni_id

print("Make the GO dictionary")
with open(go_file, 'r') as go:

    for go_line in go:
        go_line = json.loads(go_line.strip())

        if go_line["id"] not in go_dictionary:
            go_dictionary[go_line["id"]] = go_line["name"]

print("Make GO annotations dictionary")
with open(go_annotations_file, 'r') as go_annotation:

    for go_annotation_line in go_annotation:
        go_annotation_line = json.loads(go_annotation_line.strip())

        if go_annotation_line["id"] not in go_annotations_dictionary:
            go_annotations_dictionary[go_annotation_line["id"]] = []

        if go_annotation_line["go_id"] not in go_annotations_dictionary[go_annotation_line["id"]]:
            go_annotations_dictionary[go_annotation_line["id"]].append(int(go_annotation_line["go_id"].split(":")[1]))

print("Get TF's GO terms")
with open(TF_switches_file, 'r') as tf_switches:

    for tf_switches_line in tf_switches:
        tf_switches_line = tf_switches_line.strip().split('\t')
        h_tfs = tf_switches_line[2].split(",")
        non_h_tfs = tf_switches_line[3].split(",")

        for tf in h_tfs:

            if tf == "-":
                continue

            tf_uni = mapping_dictionary[tf]

            if tf_uni in go_annotations_dictionary:

                if tf_uni not in tfs_go_terms:
                    tf_uni_healthy = f"{tf_uni}_healthy"
                    tfs_go_terms[tf_uni_healthy] = []

                for go_id in go_annotations_dictionary[tf_uni]:

                    if go_dictionary[go_id] not in tfs_go_terms[tf_uni_healthy]:
                        tfs_go_terms[tf_uni_healthy].append(go_dictionary[go_id])

        for n_tf in non_h_tfs:

            if n_tf == "-":
                continue

            n_tf_uni = mapping_dictionary[n_tf]

            if n_tf_uni in go_annotations_dictionary:

                if n_tf_uni not in tfs_go_terms:
                    tf_uni_nonhealthy = f"{n_tf_uni}_nonhealthy"
                    tfs_go_terms[tf_uni_nonhealthy] = []

                for n_go_id in go_annotations_dictionary[n_tf_uni]:

                    if go_dictionary[n_go_id] not in tfs_go_terms[tf_uni_nonhealthy]:
                        tfs_go_terms[tf_uni_nonhealthy].append(go_dictionary[n_go_id])

print("Get GO term interactions")
for filename in os.listdir(input_folder):
    snp = filename.split(".")[0]

    if "statistics" in filename or "UCLeuven" in filename or "PatientSpecificNetworks" in filename:
        continue

    filepath = os.path.join(input_folder, filename)
    healthy_go_terms = []
    non_healthy_go_terms = []
    with open(filepath, 'r') as i:

        i.readline()
        for line in i:
            line = line.strip().split('\t')

            if len(line) != 2:
                continue

            healthy_go_term = line[0]
            non_healthy_go_term = line[1]

            if healthy_go_term != "-":

                if healthy_go_term not in healthy_go_terms:
                    healthy_go_terms.append(healthy_go_term)

                if healthy_go_term not in only_healthy_go_terms:
                    only_healthy_go_terms.append(healthy_go_term)

                if healthy_go_term not in all_go_terms:
                    all_go_terms.append(healthy_go_term)


            if non_healthy_go_term != "-":

                if non_healthy_go_term not in non_healthy_go_terms:
                    non_healthy_go_terms.append(non_healthy_go_term)

                if non_healthy_go_term not in only_nonhealthy_go_terms:
                    only_nonhealthy_go_terms.append(non_healthy_go_term)

                if non_healthy_go_term not in all_go_terms:
                    all_go_terms.append(non_healthy_go_term)

    for h_go_term in healthy_go_terms:

        for non_h_go_term in non_healthy_go_terms:
            go_term_interaction = f"{h_go_term}_healthy | {non_h_go_term}_nonhealthy"

            if go_term_interaction not in interactions_per_snp:
                interactions_per_snp[go_term_interaction] = []

            if snp not in interactions_per_snp[go_term_interaction]:
                interactions_per_snp[go_term_interaction].append(snp.lower())

only_h = []
only_nh = []
common_all = []
for k, v in interactions_per_snp.items():

    if len(v) < 2:
        continue

    source_i = k.split(" | ")[0].split("_")[0]
    target_i = k.split(" | ")[1].split("_")[0]

    if source_i not in only_h:
        only_h.append(source_i)

    if target_i not in only_nh:
        only_nh.append(target_i)

for xx in only_h:

    if xx in only_nh:

        if xx not in common_all:
            common_all.append(xx)

for yy in only_nh:

    if yy in only_h:

        if yy not in common_all:
            common_all.append(yy)
# print(common_all)

print("Count the appearance of GO terms in TFs (healthy and non-healthy)")
for go_term_name in all_go_terms:

    if go_term_name not in go_term_numbers_in_TFs:
        go_term_numbers_in_TFs[go_term_name] = 0

    for tf, gos in tfs_go_terms.items():

        if go_term_name in gos:
            go_term_numbers_in_TFs[go_term_name] += 1

out_nodes = []
all_stuff = []
print("Write out the interactions into the output file")
with open(output_file, 'w') as out:

    out.write("Source" + '\t' + "Target" + '\t' + "Edge.Color" + '\t' + "SNPs" + '\n')

    for int in interactions_per_snp:
        snps_all = interactions_per_snp[int]
        edge_weight = len(interactions_per_snp[int])

        if int not in all_stuff:
            all_stuff.append(int)

        # if edge_weight > 1:
        #     continue

        source_int = int.split(" | ")[0]
        target_int = int.split(" | ")[1]

        if source_int.split("_")[0] in common_all or target_int.split("_")[0] in common_all:
            continue

        out.write(source_int + '\t' + target_int + '\t' + format(edge_weight) + '\t' + ",".join(snps_all) + '\n')

        if source_int not in out_nodes:
            out_nodes.append(source_int)

        if target_int not in out_nodes:
            out_nodes.append(target_int)

print("Get the node metadata file for the network")
with open(output_node_file, 'w') as node_out:

    node_out.write("Go.Term.Name" + '\t' + "Appearance.In.TFs" + '\t' + "Type" + '\n')

    for g_term_name in out_nodes:

        only_name = g_term_name.split("_")[0]
        term_type = g_term_name.split("_")[1] 
        node_out.write(g_term_name + '\t' + format(go_term_numbers_in_TFs[only_name]) + '\t' + term_type + '\n')
