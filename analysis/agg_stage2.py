import os, re, glob, collections, statistics
OUT='results/mod_stage2'
def parse(f):
    t=open(f,errors='ignore').read()
    a=re.findall(r"'ARI':\s*([0-9.]+)",t); n=re.findall(r"'NMI':\s*([0-9.]+)",t)
    c=re.findall(r"'CA':\s*np\.float64\(([0-9.]+)\)",t) or re.findall(r"'CA':\s*([0-9.]+)",t)
    return (float(a[-1]),float(n[-1]) if n else float('nan'),float(c[-1]) if c else float('nan')) if a else None
# per (mode, dataset) mean over seeds
cell=collections.defaultdict(list)
for f in glob.glob('%s/*.log'%OUT):
    p=os.path.basename(f)[:-4].split('__')
    if len(p)!=3: continue
    mode=p[0].replace('recon_',''); ds=p[1]; r=parse(f)
    if r: cell[(mode,ds)].append(r)
modes=['both','zinb','nb','mse']
datasets=sorted({k[1] for k in cell})
mean={k:tuple(statistics.mean(x[i] for x in v) for i in range(3)) for k,v in cell.items()}
print('%-22s'%'dataset'+''.join('%9s'%m for m in modes))
for ds in datasets:
    row='%-22s'%ds[:22]
    for m in modes:
        v=mean.get((m,ds)); row+=('%9.3f'%v[0]) if v else '%9s'%'-'
    print(row)
print('\n%-10s%8s%8s%8s%10s'%('mode','ARI','NMI','ACC','fairdARI'))
for m in modes:
    aris=[mean[(m,d)][0] for d in datasets if (m,d) in mean]
    nmis=[mean[(m,d)][1] for d in datasets if (m,d) in mean]
    accs=[mean[(m,d)][2] for d in datasets if (m,d) in mean]
    # fair same-dataset delta vs both
    common=[d for d in datasets if (m,d) in mean and ('both',d) in mean]
    fd=statistics.mean(mean[(m,d)][0]-mean[('both',d)][0] for d in common) if common else 0.0
    print('%-10s%8.3f%8.3f%8.3f%+10.3f'%(m,statistics.mean(aris),statistics.mean(nmis),statistics.mean(accs),fd))
