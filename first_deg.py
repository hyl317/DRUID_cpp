# test first deg edge in DRUID graph
# I expect first-deg edge to be perfectly correct\
from collections import defaultdict

refFile = "/fs/cbsubscb09/storage/yilei/simulate/SAMAFS/safs.ped.labels"
#edgeFile = "./first_deg_edge.txt"
edgeFile = "./druid.out.1799478"
removed_inds = ['802213', 'A09222', 'A19021', 'A38060', 'A45074', 'A28102',
'840914', 'A09222', 'A27186', 'A28171', '808134']
edgeDict = defaultdict(lambda: defaultdict(lambda: "NA"))
with open(edgeFile) as edge:
    for line in edge:
        id1, id2, edge_type = line.strip().split("\t")
        edge_type = int(edge_type)
        edge_type = "PO" if edge_type == 0 else "FS"
        id1, id2 = min(id1, id2), max(id1, id2)
        edgeDict[id1][id2] = edge_type

with open(refFile) as ref:
    ref.readline() #ignore the header line
    line = ref.readline()
    while line:
        _, _, id1, id2, deg, _, rel = line.strip().split("\t")
        if deg == "Inf" or int(deg) != 1 or id1 in removed_inds or id2 in removed_inds:
            line = ref.readline()
            continue
        id1, id2 = min(id1, id2), max(id1, id2)
        if edgeDict[id1][id2] == "NA":
            print(f'missing {rel} edge between {id1} and {id2}')
        elif edgeDict[id1][id2] != rel:
            print(f'{id1} and {id2} should be {rel} but inferred as {edgeDict[id1][id2]}') 
        line = ref.readline()

