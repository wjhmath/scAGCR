import numpy as np, glob, os, re, csv
from scipy.optimize import linear_sum_assignment
PROJ="/home/liyang/BioJiaheWang/scGTAC"
def ok(cond): return "  ✓" if cond else "  ✗✗✗ 不符!"

# ---------- 工具 ----------
def parse_last_metrics(logfile):
    """从日志最后一行抓 {'CA':..,'NMI':..,'ARI':..}"""
    try: t=open(logfile).read().strip().splitlines()
    except: return None
    for line in reversed(t):
        m=re.search(r"'CA':\s*(?:np\.float64\()?([\d.]+)\)?,\s*'NMI':\s*(?:np\.float64\()?([\d.]+)\)?,\s*'ARI':\s*(?:np\.float64\()?([\d.]+)\)?", line)
        if m: return dict(ACC=float(m.group(1)),NMI=float(m.group(2)),ARI=float(m.group(3)))
    return None

# ========== 1. BENCHMARK 表1:scGTAC 七行均值 ==========
print("="*60,"\n[1] Benchmark 平均 ARI/NMI/ACC (论文 0.677/0.746/0.749)")
# scGTAC final 结果目录(按你项目:results/scagcr_final/<dataset>/...日志)
finals=sorted(glob.glob(f"{PROJ}/results/scagcr_final/*/"))
per_ds={}
for d in finals:
    ds=os.path.basename(d.rstrip("/"))
    logs=glob.glob(d+"*.log")+glob.glob(d+"*/*.log")+glob.glob(d+"*.out")
    vals=[parse_last_metrics(l) for l in logs]; vals=[v for v in vals if v]
    if vals:
        per_ds[ds]=dict(ARI=np.mean([v['ARI'] for v in vals]),
                        NMI=np.mean([v['NMI'] for v in vals]),
                        ACC=np.mean([v['ACC'] for v in vals]))
if per_ds:
    mARI=np.mean([per_ds[k]['ARI'] for k in per_ds])
    mNMI=np.mean([per_ds[k]['NMI'] for k in per_ds])
    mACC=np.mean([per_ds[k]['ACC'] for k in per_ds])
    print(f"  数据集数={len(per_ds)}  mean ARI={mARI:.3f} NMI={mNMI:.3f} ACC={mACC:.3f}")
    print("  vs 论文 0.677/0.746/0.749:",ok(abs(mARI-0.677)<.005 and abs(mNMI-0.746)<.005 and abs(mACC-0.749)<.005))
else:
    print("  [未找到 results/scagcr_final/*/ 日志 —— 请把你 scGTAC 最终结果目录名告诉我,我改路径]")

# ========== 2. Muraro 案例:ARI 0.913 + 内分泌召回(表2) ==========
print("="*60,"\n[2] Muraro 案例 (论文 ARI 0.913; α0.97 β0.95 δ0.97 PP1.00)")
emb=f"{PROJ}/paper_figures/Figure3_umap/emb_muraro_scGTAC.npz"
if os.path.exists(emb):
    z=np.load(emb,allow_pickle=True)
    y=list(map(str,np.asarray(z['y']).ravel())); p=list(map(str,np.asarray(z['pred']).ravel()))
    from sklearn.metrics import adjusted_rand_score
    print(f"  scGTAC Muraro ARI={adjusted_rand_score(y,p):.3f}  vs 0.913:",ok(abs(adjusted_rand_score(y,p)-0.913)<.005))
    def recall_by_type(true,pred):
        types=sorted(set(true)); clus=sorted(set(pred)); C=np.zeros((len(types),len(clus)))
        for a,b in zip(true,pred): C[types.index(a),clus.index(b)]+=1
        r,cc=linear_sum_assignment(-C); rec={}
        for i,j in zip(r,cc): rec[types[i]]=C[i,j]/C[i].sum() if C[i].sum()>0 else 0
        return rec
    rec=recall_by_type(y,p)
    want={"pancreatic A cell":0.97,"type B pancreatic cell":0.95,"pancreatic D cell":0.97,"pancreatic PP cell":1.00}
    for t,w in want.items():
        got=rec.get(t,float('nan')); print(f"  {t:26s} recall={got:.2f} vs {w}",ok(abs(got-w)<.02))
else:
    print("  [未找到 emb_muraro_scGTAC.npz]")

# ========== 3. 表2 baseline 内分泌召回 ==========
print("="*60,"\n[3] 表2 baseline 召回 (Seurat α0.67/β0.61; scGNN δ0.00; scVI β0.50; DEC α0.78; scDSC PP0.00)")
def recall_named(true_codes,pred,ynames):
    c2n={}; 
    for c,n in zip(true_codes,ynames): c2n[c]=n
    true=[c2n[c] for c in true_codes]
    types=sorted(set(true)); clus=sorted(set(pred)); C=np.zeros((len(types),len(clus)))
    for a,b in zip(true,pred): C[types.index(a),clus.index(b)]+=1
    r,cc=linear_sum_assignment(-C); rec={}
    for i,j in zip(r,cc): rec[types[i]]=C[i,j]/C[i].sum() if C[i].sum()>0 else 0
    return rec
if os.path.exists(emb):
    ynames=list(map(str,np.asarray(z['y']).ravel())); N=len(ynames)
    checks={"seurat":("pancreatic A cell",0.67),"scgnn":("pancreatic D cell",0.00),
            "scvi":("type B pancreatic cell",0.50),"dec":("pancreatic A cell",0.78),
            "scdsc":("pancreatic PP cell",0.00)}
    import pandas as pd
    for m,(typ,w) in checks.items():
        fs=glob.glob(f"{PROJ}/results/baselines/{m}/muraro_pancreas_seed1_pred.csv")
        if not fs: print(f"  [{m} 缺 pred.csv]"); continue
        df=pd.read_csv(fs[0])
        if len(df)!=N: print(f"  [{m} 行数 {len(df)}≠{N}]"); continue
        rec=recall_named(list(map(str,df['true'])),list(map(str,df['pred'])),ynames)
        got=rec.get(typ,float('nan')); print(f"  {m:8s} {typ:24s}={got:.2f} vs {w}",ok(abs(got-w)<.03))

# ========== 4. 消融 15-set ==========
print("="*60,"\n[4] 消融 15-set Δ(full−variant) (论文 CL ARI0.115/ACC0.095; graph .022/.019; edge .015/.015; ZINB .005/.001)")
EXCLUDE={"10X_PBMC","zheng68k"}
def variant_means(v):
    ds=sorted(d for d in glob.glob(f"{PROJ}/results/ablation/{v}/*/") if os.path.basename(d.rstrip('/')) not in EXCLUDE)
    A,Nn,Cc=[],[],[]
    for d in ds:
        logs=glob.glob(d+"*.log")+glob.glob(d+"*.out")
        vals=[parse_last_metrics(l) for l in logs]; vals=[x for x in vals if x]
        if vals: A.append(np.mean([x['ARI'] for x in vals]));Nn.append(np.mean([x['NMI'] for x in vals]));Cc.append(np.mean([x['ACC'] for x in vals]))
    return (np.mean(A),np.mean(Nn),np.mean(Cc),len(A)) if A else (None,)*4
full=variant_means("full")
if full[0] is not None:
    print(f"  full: ARI={full[0]:.3f} NMI={full[1]:.3f} ACC={full[2]:.3f} (n={full[3]})")
    want={"wo_cl":(0.115,0.095),"wo_graph":(0.022,0.019),"wo_edge":(0.015,0.015),"wo_recon":(0.005,0.001)}
    for v,(dA,dC) in want.items():
        vm=variant_means(v)
        if vm[0] is None: print(f"  [{v} 缺]"); continue
        gA,gC=full[0]-vm[0], full[2]-vm[2]
        print(f"  {v:9s} ΔARI={gA:+.3f}(vs {dA}) ΔACC={gC:+.3f}(vs {dC})",ok(abs(gA-dA)<.01 and abs(gC-dC)<.01))
else:
    print("  [未找到 results/ablation/full/*/ —— 路径可能不同]")

# ========== 5. 鲁棒/可扩展 ==========
print("="*60,"\n[5] 鲁棒(Baron 10% >0.89; UUO full≈0.61) + 可扩展(68k≈288s)")
ds_csv=f"{PROJ}/results/robustness/downsample.csv"
if os.path.exists(ds_csv):
    rows=[r for r in csv.reader(open(ds_csv))]
    b10=[float(r[3]) for r in rows if len(r)>3 and r[0].lower()=="baron" and r[1] in("0.1","0.10")]
    if b10: print(f"  Baron 10% ARI 均值={np.mean(b10):.3f}",ok(np.mean(b10)>0.89))
tcsv=f"{PROJ}/results/scalability/timing.csv"
if os.path.exists(tcsv):
    for r in csv.reader(open(tcsv)):
        if r and "68k" in r[0].lower().replace("pbmc",""):
            print(f"  68k 计时={r[-1]}s",ok(200<float(r[-1])<400))
print("="*60,"\n核查完成。把整段输出贴给我。")
