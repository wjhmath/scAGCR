import os, re, glob, collections, statistics
OUT='results/mod_featv2'
def parse(f):
    t=open(f,errors='ignore').read()
    a=re.findall(r"'ARI':\s*([0-9.]+)",t); n=re.findall(r"'NMI':\s*([0-9.]+)",t)
    c=re.findall(r"'CA':\s*np\.float64\(([0-9.]+)\)",t) or re.findall(r"'CA':\s*([0-9.]+)",t)
    return (float(a[-1]),float(n[-1]) if n else float('nan'),float(c[-1]) if c else float('nan')) if a else None
agg=collections.defaultdict(list)
for f in glob.glob('%s/*.log'%OUT):
    p=os.path.basename(f)[:-4].split('__')
    if len(p)!=3: continue
    r=parse(f)
    if r: agg[p[0]].append(r)
rows=[(t,len(v),statistics.mean(x[0] for x in v),statistics.mean(x[1] for x in v),statistics.mean(x[2] for x in v)) for t,v in agg.items()]
b=dict((r[0],r[2]) for r in rows).get('feat_none')
print('%-14s%4s%8s%8s%8s%9s'%('strategy','n','ARI','NMI','ACC','dARI_none'))
for t,nn,a,n,c in sorted(rows,key=lambda x:-x[2]):
    print('%-14s%4d%8.3f%8.3f%8.3f%+9.3f'%(t,nn,a,n,c,(a-b) if b is not None else 0))
