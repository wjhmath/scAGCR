import os, re, glob, collections, statistics
EXCLUDE={'10X_PBMC','zheng68k'}
def parse(f):
    a=re.findall(r"'ARI':\s*([0-9.]+)",open(f,errors='ignore').read())
    return float(a[-1]) if a else None
def load(d,tagfilter=None,rename=None):
    cell=collections.defaultdict(list)
    for f in glob.glob('%s/*.log'%d):
        p=os.path.basename(f)[:-4].split('__')
        if len(p)!=3 or p[1] in EXCLUDE: continue
        if tagfilter and p[0]!=tagfilter: continue
        r=parse(f)
        if r is not None: cell[((rename or p[0]),p[1],p[2])].append(r)
    return cell
none=load('results/mod_abl2','wo_feat','none')      # none(zinb) 3 seeds 已有
gnz=load('results/mod_featconf')
all={**none,**gnz}
# per (tag,dataset) mean over seeds
m=collections.defaultdict(dict)
for (tag,ds,sd),v in all.items(): m[tag][ds]=m[tag].get(ds,[])+v
mean={(tag,ds):statistics.mean(vs) for tag,d in m.items() for ds,vs in d.items()}
ds_none=set(d for (t,d) in mean if t=='none')
print('%-10s%8s%9s%9s'%('strategy','ARI','dvs_none','wins/15'))
nmean=statistics.mean(mean[('none',d)] for d in ds_none)
print('%-10s%8.3f%9s%9s'%('none',nmean,'-','-'))
for tag in ['gnz005','gnz01','gnz015']:
    ds=[d for d in ds_none if (tag,d) in mean]
    if not ds: continue
    a=statistics.mean(mean[(tag,d)] for d in ds)
    dl=statistics.mean(mean[(tag,d)]-mean[('none',d)] for d in ds)
    wins=sum(1 for d in ds if mean[(tag,d)]-mean[('none',d)]>0.002)
    print('%-10s%8.3f%+9.3f%6d/%d'%(tag,a,dl,wins,len(ds)))
