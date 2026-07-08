import os, re, glob, collections, statistics
def parse(f):
    t=open(f,errors='ignore').read()
    a=re.findall(r"'ARI':\s*([0-9.]+)",t); n=re.findall(r"'NMI':\s*([0-9.]+)",t)
    c=re.findall(r"'CA':\s*np\.float64\(([0-9.]+)\)",t) or re.findall(r"'CA':\s*([0-9.]+)",t)
    return (float(a[-1]),float(n[-1]) if n else float('nan'),float(c[-1]) if c else float('nan')) if a else None
cell=collections.defaultdict(list)
# lambda 0.2 from stage2 zinb
for f in glob.glob('results/mod_stage2/recon_zinb__*.log'):
    p=os.path.basename(f)[:-4].split('__'); r=parse(f)
    if len(p)==3 and r: cell[('0.2',p[1])].append(r)
for f in glob.glob('results/mod_lambda/*.log'):
    p=os.path.basename(f)[:-4].split('__')
    if len(p)!=3: continue
    lam={'l050':'0.5','l100':'1.0','l200':'2.0'}.get(p[0].replace('zinb_',''))
    r=parse(f)
    if lam and r: cell[(lam,p[1])].append(r)
lams=['0.2','0.5','1.0','2.0']; datasets=sorted({k[1] for k in cell})
mean={k:tuple(statistics.mean(x[i] for x in v) for i in range(3)) for k,v in cell.items()}
print('%-10s%8s%8s%8s%10s'%('lambda','ARI','NMI','ACC','fairdARI'))
for L in lams:
    aris=[mean[(L,d)][0] for d in datasets if (L,d) in mean]
    nmis=[mean[(L,d)][1] for d in datasets if (L,d) in mean]
    accs=[mean[(L,d)][2] for d in datasets if (L,d) in mean]
    common=[d for d in datasets if (L,d) in mean and ('0.2',d) in mean]
    fd=statistics.mean(mean[(L,d)][0]-mean[('0.2',d)][0] for d in common) if common else 0.0
    if aris: print('%-10s%8.3f%8.3f%8.3f%+10.3f'%(L,statistics.mean(aris),statistics.mean(nmis),statistics.mean(accs),fd))
