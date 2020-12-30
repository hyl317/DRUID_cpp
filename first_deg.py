# test first deg edge in DRUID graph
# I expect first-deg edge to be perfectly correct\
from collections import defaultdict

refFile = "/fs/cbsubscb09/storage/yilei/simulate/SAMAFS/safs.ped.labels"
edgeFile = "./safs.DRUID"
removed_inds = ['802213', 'A09222', 'A19021', 'A38060', 'A45074', 'A28102',
'840914', 'A09222', 'A27186', 'A28171', '808134']
edgeDict = defaultdict(lambda: defaultdict(lambda: "NA"))

edgeMap = {0: "PO", 1:"FS", 2:"GP", 3:"AV"}
with open(edgeFile) as edge:
    for line in edge:
        id1, id2, edge_type = line.strip().split("\t")
        if edge_type not in ['PC', 'AV', 'FS']:
            continue
        id1, id2 = min(id1, id2), max(id1, id2)
        if edge_type == 'PC':
            edge_type = 'PO'
        edgeDict[id1][id2] = edge_type

with open(refFile) as ref:
    ref.readline() #ignore the header line
    line = ref.readline()
    while line:
        _, _, id1, id2, deg, _, rel = line.strip().split("\t")
        if deg == "Inf" or int(deg) > 2 or id1 in removed_inds or id2 in removed_inds:
            line = ref.readline()
            continue
        id1, id2 = min(id1, id2), max(id1, id2)
        #if edgeDict[id1][id2] == "NA":
        #    print(f'missing {rel} edge between {id1} and {id2}')
        if edgeDict[id1][id2] != rel and edgeDict[id1][id2] != "NA":
            print(f'{id1} and {id2} should be {rel} but inferred as {edgeDict[id1][id2]}') 
        line = ref.readline()

