import os, re, glob, collections, statistics
OUT='results/mod_innov2'
def parse(f):
    a=re.findall(r"'ARI':\s*([0-9.]+)",open(f,errors='ignore').read())
    return float(a[-1]) if a else None
agg=collections.defaultdict(list)
for f in glob.glob('%s/*.log'%OUT):
    p=os.path.basename(f)[:-4].split('__')
    if len(p)!=3: continue
    r=parse(f)
    if r is not None: agg[p[0]].append(r)
rows=[(t,len(v),statistics.mean(v)) for t,v in agg.items()]
b=dict((t,a) for t,n,a in rows).get('base')
print('!!! DEV ONLY — 越线者必须过全15关(>=+0.010 且 >=10/15)!!!')
print('%-10s%4s%8s%9s'%('config','n','ARI','dvs_base'))
for t,n,a in sorted(rows,key=lambda x:-x[2]):
    print('%-10s%4d%8.3f%+9.3f'%(t,n,a,(a-b) if b is not None else 0))
