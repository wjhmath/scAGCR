import os, re, glob, collections, statistics
OUT='results/mod_abl2'; EXCLUDE={'10X_PBMC','zheng68k'}
def parse(f):
    a=re.findall(r"'ARI':\s*([0-9.]+)",open(f,errors='ignore').read())
    return float(a[-1]) if a else None
cell=collections.defaultdict(list)
for f in glob.glob('%s/*.log'%OUT):
    p=os.path.basename(f)[:-4].split('__')
    if len(p)!=3 or p[1] in EXCLUDE: continue
    r=parse(f)
    if r is not None: cell[(p[0],p[1])].append(r)
order=['full','wo_cl','wo_recon','wo_feat','wo_edge','wo_graph']
ds=sorted({k[1] for k in cell}); mean={k:statistics.mean(v) for k,v in cell.items()}
print('=== Δ vs full(负=该组件有用)===')
print('%-20s'%'dataset'+''.join('%9s'%o.replace('wo_','-') for o in order[1:]))
for d in ds:
    row='%-20s'%d[:20]; fu=mean.get(('full',d))
    for o in order[1:]:
        row+=('%+9.3f'%(mean[(o,d)]-fu)) if ((o,d) in mean and fu is not None) else '%9s'%'-'
    print(row)
