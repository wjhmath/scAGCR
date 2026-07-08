import os, re, glob, collections, statistics
EXCLUDE={'10X_PBMC','zheng68k'}   # 论文的两个 scalability 附加集
def parse(f):
    t=open(f,errors='ignore').read()
    a=re.findall(r"'ARI':\s*([0-9.]+)",t); n=re.findall(r"'NMI':\s*([0-9.]+)",t)
    c=re.findall(r"'CA':\s*np\.float64\(([0-9.]+)\)",t) or re.findall(r"'CA':\s*([0-9.]+)",t)
    return (float(a[-1]),float(n[-1]) if n else float('nan'),float(c[-1]) if c else float('nan')) if a else None
def load(d, tagfix=lambda s:s):
    cell=collections.defaultdict(list)
    for f in glob.glob('%s/*.log'%d):
        p=os.path.basename(f)[:-4].split('__')
        if len(p)!=3 or p[1] in EXCLUDE: continue
        r=parse(f)
        if r: cell[(tagfix(p[0]),p[1])].append(r)
    return cell
def table(cell, order, base, title):
    mean={k:tuple(statistics.mean(x[i] for x in v) for i in range(3)) for k,v in cell.items()}
    dsall=sorted({k[1] for k in cell})
    print('\n===== %s (15 数据集) ====='%title)
    print('%-10s %7s %7s %7s | %7s %7s %7s'%('variant','ARI','NMI','ACC','dARI','dNMI','dACC'))
    for t in order:
        ds=[d for d in dsall if (t,d) in mean]
        if not ds: continue
        raw=tuple(statistics.mean(mean[(t,d)][i] for d in ds) for i in range(3))
        common=[d for d in dsall if (t,d) in mean and (base,d) in mean]
        dl=tuple(statistics.mean(mean[(t,d)][i]-mean[(base,d)][i] for d in common) for i in range(3)) if common else (0,0,0)
        print('%-10s %7.3f %7.3f %7.3f | %+7.3f %+7.3f %+7.3f'%(t,raw[0],raw[1],raw[2],dl[0],dl[1],dl[2]))
# 主结果:重建模式对比
table(load('results/mod_stage2', lambda s:s.replace('recon_','')), ['both','zinb','nb','mse'], 'both', '主结果 重建对比 (Δ vs both=原版)')
# 消融:新方法(zinb-only)
table(load('results/mod_abl2'), ['full','wo_cl','wo_recon','wo_feat','wo_edge','wo_graph'], 'full', '消融 (Δ vs full=新方法)')
