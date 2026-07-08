import numpy as np, glob, os, json, re
PROJ="/home/liyang/BioJiaheWang/scGTAC"
LINEUP=["muraro_pancreas","baron","GSE103354","GSE103322","Tonsil","GSE150580_Mammary",
        "GSE119531","GSE194122_PBMC_Bench_1","multiome","GSE159115_ccRCC",
        "GSE194122_PBMC_Test","Goolam","Crohn","68kPBMC","GSE123516_labeled"]
BASE=["seurat","scgnn","scdeepcluster","scvi","dec","scdsc"]
DISP={"scGTAC":"scGTAC(ours)","seurat":"Seurat","scgnn":"scGNN","scdeepcluster":"scDeepCluster",
      "scvi":"scVI","dec":"DEC","scdsc":"scDSC"}
def scg(ds):
    v={"ARI":[],"NMI":[],"CA":[]}
    for f in glob.glob(f"{PROJ}/results/scagcr_final/{ds}/run_seed*.log"):
        t=open(f,errors="ignore").read().strip().splitlines()
        if not t: continue
        for k in v:
            m=re.search(rf"'{k}':\s*(?:np\.float64\()?([0-9.]+)",t[-1])
            if m: v[k].append(float(m.group(1)))
    return (np.mean(v["ARI"]) if v["ARI"] else np.nan, np.mean(v["NMI"]) if v["NMI"] else np.nan, np.mean(v["CA"]) if v["CA"] else np.nan)
def base(m,ds):
    A,N,C=[],[],[]
    for f in glob.glob(f"{PROJ}/results/baselines/{m}/{ds}_seed*_metrics.json")+glob.glob(f"{PROJ}/results/baselines_chuli/{m}/{ds}_seed*_metrics.json"):
        try:
            j=json.load(open(f)); A.append(j["ARI"]); N.append(j["NMI"]); C.append(j["ACC"])
        except Exception: pass
    return (np.mean(A) if A else np.nan, np.mean(N) if N else np.nan, np.mean(C) if C else np.nan)
order=["scGTAC"]+BASE
M={"scGTAC":{ds:scg(ds) for ds in LINEUP}}
for m in BASE: M[m]={ds:base(m,ds) for ds in LINEUP}
print("="*64);print("CHECK4  方法均值表 (15 数据集均值;对照正文表)");print("="*64)
print("  方法             ARI    NMI    ACC   缺数据集")
for m in order:
    arr=np.array([M[m][ds] for ds in LINEUP]); mean=np.nanmean(arr,axis=0)
    nmiss=int(np.isnan(arr[:,0]).sum())
    print("  %-15s %.3f  %.3f  %.3f   %s"%(DISP[m],mean[0],mean[1],mean[2], (f"{nmiss}个" if nmiss else "")))
print("  论文表: Seurat 0.499/0.680/0.610 ; scGTAC 0.677/0.746/0.749")
print("="*64);print("CHECK4b  scGTAC 第一名计数 (15×3=45)");print("="*64)
cnt=0; tot=0; lose=[]
for ds in LINEUP:
    for mi,mn in enumerate(["ARI","NMI","ACC"]):
        vals={m:M[m][ds][mi] for m in order if not np.isnan(M[m][ds][mi])}
        if "scGTAC" not in vals: continue
        tot+=1
        if vals["scGTAC"]>=max(vals.values())-1e-9: cnt+=1
        else:
            win=max(vals,key=vals.get); lose.append(f"{ds}/{mn}(scGTAC {vals['scGTAC']:.3f} < {DISP[win]} {vals[win]:.3f})")
print(f"  scGTAC 第一名: {cnt}/{tot}   (论文声明 39/45)")
if lose:
    print("  非第一名的格子:")
    for l in lose: print("    -",l)
