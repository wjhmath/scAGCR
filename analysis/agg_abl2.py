import os, re, glob, collections, statistics
OUT='results/mod_abl2'
def parse(f):
    t=open(f,errors='ignore').read()
    a=re.findall(r"'ARI':\s*([0-9.]+)",t); n=re.findall(r"'NMI':\s*([0-9.]+)",t)
    c=re.findall(r"'CA':\s*np\.float64\(([0-9.]+)\)",t) or re.findall(r"'CA':\s*([0-9.]+)",t)
    return (float(a[-1]),float(n[-1]) if n else float('nan'),float(c[-1]) if c else float('nan')) if a else None
cell=collections.defaultdict(list)
for f in glob.glob('%s/*.log'%OUT):
    p=os.path.basename(f)[:-4].split('__')
    if len(p)!=3: continue
    r=parse(f)
    if r: cell[(p[0],p[1])].append(r)
order=['full','wo_cl','wo_recon','wo_feat','wo_edge','wo_graph']
datasets=sorted({k[1] for k in cell})
mean={k:tuple(statistics.mean(x[i] for x in v) for i in range(3)) for k,v in cell.items()}
print('%-10s%8s%8s%8s%10s'%('variant','ARI','NMI','ACC','dARI_vs_full'))
for t in order:
    aris=[mean[(t,d)][0] for d in datasets if (t,d) in mean]
    nmis=[mean[(t,d)][1] for d in datasets if (t,d) in mean]
    accs=[mean[(t,d)][2] for d in datasets if (t,d) in mean]
    common=[d for d in datasets if (t,d) in mean and ('full',d) in mean]
    fd=statistics.mean(mean[(t,d)][0]-mean[('full',d)][0] for d in common) if common else 0.0
    if aris: print('%-10s%8.3f%8.3f%8.3f%+10.3f'%(t,statistics.mean(aris),statistics.mean(nmis),statistics.mean(accs),fd))
