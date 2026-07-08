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
print('%-10s | %-22s | %-22s'%('','        raw mean','     fair Δ vs full'))
print('%-10s %7s %7s %7s | %7s %7s %7s'%('variant','ARI','NMI','ACC','dARI','dNMI','dACC'))
for t in order:
    ds=[d for d in datasets if (t,d) in mean]
    raw=tuple(statistics.mean(mean[(t,d)][i] for d in ds) for i in range(3))
    common=[d for d in datasets if (t,d) in mean and ('full',d) in mean]
    dl=tuple(statistics.mean(mean[(t,d)][i]-mean[('full',d)][i] for d in common) for i in range(3)) if common else (0,0,0)
    print('%-10s %7.3f %7.3f %7.3f | %+7.3f %+7.3f %+7.3f'%(t,raw[0],raw[1],raw[2],dl[0],dl[1],dl[2]))
