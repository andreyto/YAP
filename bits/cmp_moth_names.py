if False:
    g_sets = []
    for g in gr:
        print g
        g_seq = open(g,"r").readlines()
        g_s = [ x.split()[0].strip() for x in g_seq ]
        g_sets.append((g_seq[0].strip().split()[1],set(g_s)))

    n_sets = []
    for n in nm:
        n_s = open(n,"r").readlines()
        n_s = [ l.strip().split()[1].split(",") for l in n_s ]
        n_set = set()
        for n_s_e in n_s:
            n_set |= set(n_s_e)
        print "{} len {}".format(n,len(n_set))
        n_sets.append((n,n_set))


for n in n_sets:
    print n[0],":",
    for g in g_sets:
        if g[1] & n[1]:
            print g[0]
            break
    else:
        print ""

if False:

    for g in g_sets:
        print g[0],":",
        for n in n_sets:
            if g[1] & n[1]:
                print n[0]
                break
        else:
            print ""


