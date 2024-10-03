input_file = "pipeline_results_11_04_2024/UC/UC_patients_with_clusters_13_04_2024.tsv"
output_file = "pipeline_results_11_04_2024/UC/UC_number_of_goterms_per_patients_13_04_2024.tsv"

patients_dictionary = {}
all_go_terms = {}

with open(input_file, 'r') as i:

    i.readline()

    for line in i:
        line = line.strip().split('\t')

        healthy_goterm = line[0].split("|")[0]
        non_healthy_goterm = line[0].split("|")[1]

        for x in range(1, len(line)):
            patient_id = f"Patient_{x}"

            if patient_id not in patients_dictionary:
                patients_dictionary[patient_id] = {"healthy": [], "non_healthy": []}

            if line[x] == "0":
                continue

            else:

                if healthy_goterm not in patients_dictionary[patient_id]["healthy"]:
                    patients_dictionary[patient_id]["healthy"].append(healthy_goterm)

                if healthy_goterm not in all_go_terms:
                    all_go_terms[healthy_goterm] = {"healthy": 0, "non_healthy": 0}

                if non_healthy_goterm not in patients_dictionary[patient_id]["non_healthy"]:
                    patients_dictionary[patient_id]["non_healthy"].append(non_healthy_goterm)

                if non_healthy_goterm not in all_go_terms:
                    all_go_terms[non_healthy_goterm] = {"healthy": 0, "non_healthy": 0}

for goterm in all_go_terms:

    for k, v in patients_dictionary.items():

        if goterm in v["healthy"]:
            all_go_terms[goterm]["healthy"] += 1

        if goterm in v["non_healthy"]:
            all_go_terms[goterm]["non_healthy"] += 1

with open(output_file, 'w') as out:
    out.write("GOTermName" + '\t' + "NumberOfPatients(Healthy)" + '\t' + "NumberOfPatients(Non-healthy)" + '\n')

    for x, y in all_go_terms.items():

        out.write(x + '\t' + format(y["healthy"]) + '\t' + format(y["non_healthy"]) + '\n')
