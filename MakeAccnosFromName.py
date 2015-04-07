import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("-n","--name", help="name file")
parser.add_argument("-s","--min-cluster-size", type=int, help="cluster size cutoff")
args = parser.parse_args()

name = args.name.strip()
assert len(name) > 0
accnos = os.path.basename(name).rsplit(".name",1)[0]+".accnos"
with open(name,"r") as inp, open(accnos,"w") as out:
    for line in inp:
        line = line.strip()
        if line:
            rep, clust = line.split("\t")
            clust = [ y for y in ( x.strip() for x in clust.split(",") ) if y ]
            if len(clust) < args.min_cluster_size:
                ## add all sequences in the cluster to be sure
                ## clust shold include rep already but just in case
                out.write("\n".join(list(set(clust)|set([rep])))+"\n")


