import os, re, glob, collections, statistics
EXCLUDE={'10X_PBMC','zheng68k'}
def parse(f):
    a=re.findall(r"'ARI':\s*([0-9.]+)",open(f,errors='ignore').read())
    return float(a[-1]) if a else None
def load(d,tagfilter=None,rename=None):
    c=collections.defaultdict(list)
    for f in glob.glob('%s/*.log'%d):
        p=os.path.basename(f)[:-4].split('__')
        if len(p)!=3 or p[1] in EXCLUDE: continue
        if tagfilter and p[0]!=tagfilter: continue
        r=parse(f)
        if r is not None: c[((rename or p[0]),p[1])].append(r)
    return c
allc={}
allc.update(load('results/mod_abl2','wo_feat','none'))
allc.update(load('results/mod_featconf2'))
m=collections.defaultdict(dict)
for (t,ds),v in allc.items(): m[t][ds]=statistics.mean(v)
dsn=set(m['none'])
print('%-10s%8s%9s%9s   判定'%('cand','ARI','dvs_none','wins/15'))
print('%-10s%8.3f%9s%9s'%('none',statistics.mean(m['none'].values()),'-','-'))
for t in ['mask015','mask005','gauss002']:
    if t not in m: continue
    ds=[d for d in dsn if d in m[t]]
    a=statistics.mean(m[t][d] for d in ds)
    dl=statistics.mean(m[t][d]-m['none'][d] for d in ds)
    wins=sum(1 for d in ds if m[t][d]-m['none'][d]>0.002)
    ok='✓真增强' if (dl>=0.010 and wins>=10) else '✗噪声'
    print('%-10s%8.3f%+9.3f%6d/%d   %s'%(t,a,dl,wins,len(ds),ok))
