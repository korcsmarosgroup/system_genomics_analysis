import json
import os


mapping_file = "files/human_uniprot_mapping_version_2022_01.json"
input_folder = "pipeline_results_11_04_2024/CD/cytoscape_networks"
tf_switches = "pipeline_results_11_04_2024/CD/TF_switches_results_17_11_2021_CDLeuven_processed.tsv"
mapping_dictionary = {}
tfnames_dictionary = {}

with open(mapping_file, 'r') as mapping:

    mapping.readline()

    for mapping_line in mapping:
        mapping_line = json.loads(mapping_line.strip())

        if mapping_line["from_id_type"] == "genename_primary" and mapping_line["to_id_type"] == "uniprotac":
            gene_name = mapping_line["from_id"].upper()
            uniprotac = mapping_line["to_id"].upper()

            if uniprotac not in mapping_dictionary:
                mapping_dictionary[uniprotac] = gene_name

with open(tf_switches, 'r') as tf_switches:

    for tf_line in tf_switches:
        tf_line = tf_line.strip().split('\t')

        if tf_line[0] not in tfnames_dictionary:
            tfnames_dictionary[tf_line[0]] = {"healthy": [], "disease": []}

        for source_tf in tf_line[2].split(","):

            if source_tf not in tfnames_dictionary[tf_line[0]]["healthy"]:
                tfnames_dictionary[tf_line[0]]["healthy"].append(source_tf)

        for target_tf in tf_line[3].split(","):

            if target_tf not in tfnames_dictionary[tf_line[0]]["disease"]:
                tfnames_dictionary[tf_line[0]]["disease"].append(target_tf)

for filename in os.listdir(input_folder):

    if filename == "final_networks":
        continue

    if "CD" in filename:
        print(filename)

        actual_file = os.path.join("pipeline_results_11_04_2024/CD/cytoscape_networks", filename)
        actual_filename = filename.split(".")[0]
        new_filename = f"{actual_filename}_final.tsv"
        output_file = os.path.join("pipeline_results_11_04_2024/CD/cytoscape_networks/final_networks", new_filename)

        with open(actual_file, 'r') as inp, open(output_file, 'w') as out:
            out.write("HealthyGOTerm\tDiseaseGOTerm\tHealthyGOTermTFs\tDiseaseGOTermTFs\tRelevantSNPs\n")
            inp.readline()

            for line in inp:
                line = line.strip().split('\t')
                source = line[0]
                target = line[1]
                source_tfs = line[2].split(",")
                target_tfs = line[3].split(",")
                snps = line[4].split(",")
                
                healthy_tfs = []
                disease_tfs = []

                for s in snps:
                    s = s.upper()

                    for sourcetf in source_tfs:

                        if sourcetf in mapping_dictionary:

                            if mapping_dictionary[sourcetf] in tfnames_dictionary[s]["healthy"]:

                                if mapping_dictionary[sourcetf] not in healthy_tfs:
                                    healthy_tfs.append(mapping_dictionary[sourcetf])

                    for targettf in target_tfs:

                        if targettf in mapping_dictionary:

                            if mapping_dictionary[targettf] in tfnames_dictionary[s]["disease"]:

                                if mapping_dictionary[targettf] not in disease_tfs:
                                    disease_tfs.append(mapping_dictionary[targettf])

                healthy_tfs_out = ",".join(healthy_tfs)
                disease_tfs_out = ",".join(disease_tfs)
                snps_out = ",".join(snps)

                out.write(f"{source}\t{target}\t{healthy_tfs_out}\t{disease_tfs_out}\t{snps_out}\n")
