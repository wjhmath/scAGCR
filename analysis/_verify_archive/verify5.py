import numpy as np, glob, pandas as pd
from scipy.optimize import linear_sum_assignment
PROJ="/home/liyang/BioJiaheWang/scGTAC"
z=np.load(f"{PROJ}/paper_figures/Figure3_umap/emb_muraro_scGTAC.npz",allow_pickle=True)
ytrue=list(map(str,np.asarray(z['y']).ravel())); N=len(ytrue)
def recall_by_type(true,pred):
    types=sorted(set(true)); clus=sorted(set(pred))
    C=np.zeros((len(types),len(clus)))
    for a,b in zip(true,pred): C[types.index(a),clus.index(b)]+=1
    r,cc=linear_sum_assignment(-C); rec={t:0.0 for t in types}
    for i,j in zip(r,cc):
        if C[i].sum()>0: rec[types[i]]=C[i,j]/C[i].sum()
    return rec
DISP={"scGTAC":"scGTAC","seurat":"Seurat","scgnn":"scGNN","scdeepcluster":"scDeepC",
      "scvi":"scVI","dec":"DEC","scdsc":"scDSC"}
res={"scGTAC":recall_by_type(ytrue, list(map(str,np.asarray(z['pred']).ravel())))}
for m in ["seurat","scgnn","scdeepcluster","scvi","dec","scdsc"]:
    fs=glob.glob(f"{PROJ}/results/baselines/{m}/muraro_pancreas_seed1_pred.csv")+\
       glob.glob(f"{PROJ}/results/baselines_chuli/{m}/muraro_pancreas_seed1_pred.csv")
    if not fs: print(f"  [缺] {m}"); continue
    df=pd.read_csv(fs[0])
    pcol="pred" if "pred" in df.columns else df.columns[0]
    tcol="true" if "true" in df.columns else (df.columns[1] if len(df.columns)>1 else None)
    pred=list(map(str,df[pcol]))
    if tcol is not None:
        traw=list(map(str,df[tcol]))
        if len(df)==N:
            c2n={}
            for c,n in zip(traw,ytrue): c2n[c]=n
            true=[c2n[c] for c in traw]
        else: true=traw
    elif len(df)==N:
        true=ytrue
    else:
        print(f"  [跳过] {m}: 列={list(df.columns)} 行={len(df)} vs N={N}"); continue
    res[m]=recall_by_type(true,pred)
types=sorted(set(ytrue)); order=[m for m in ["scGTAC","seurat","scgnn","scdeepcluster","scvi","dec","scdsc"] if m in res]
print("逐亚型 recall (对照 Fig3b):")
print("  %-26s"%"cell type"+"".join("%-9s"%DISP[m] for m in order))
for t in types:
    print("  %-26s"%t+"".join("%-9.2f"%res[m].get(t,np.nan) for m in order))
