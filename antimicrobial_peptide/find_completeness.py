query = "merged_results_ids.txt"
reference = "concat.okay.fa.completeness.list"

reference_dict = dict()

with open(reference, "r") as f:
    for line in f:
        line = line.split("\t")
        # print(line)
        reference_dict[line[0]] = line[1].rstrip()

# print(reference_dict)

with open(query, "r") as f:
    for line in f:
        output = line.rstrip() + "\t" + reference_dict[line.rstrip()]
        print(output)