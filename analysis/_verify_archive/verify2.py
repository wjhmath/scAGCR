import numpy as np, glob, re, os
PROJ="/home/liyang/BioJiaheWang/scGTAC"
LINEUP=["muraro_pancreas","baron","GSE103354","GSE103322","Tonsil","GSE150580_Mammary",
        "GSE119531","GSE194122_PBMC_Bench_1","multiome","GSE159115_ccRCC",
        "GSE194122_PBMC_Test","Goolam","Crohn","68kPBMC","GSE123516_labeled"]
EX={"10X_PBMC","zheng68k"}
def rd(d):
    v={"ARI":[],"NMI":[],"CA":[]}
    for f in glob.glob(d+"/run_seed*.log"):
        t=open(f,errors="ignore").read().strip().splitlines()
        if not t: continue
        for k in v:
            m=re.search(rf"'{k}':\s*(?:np\.float64\()?([0-9.]+)",t[-1])
            if m: v[k].append(float(m.group(1)))
    return {k:(np.mean(x) if x else np.nan) for k,x in v.items()}
print("="*62);print("CHECK1  benchmark 15-set 均值  (论文: 0.677 / 0.746 / 0.749)");print("="*62)
per={d:rd(f"{PROJ}/results/scagcr_final/{d}") for d in LINEUP}
mi=[d for d in LINEUP if np.isnan(per[d]["ARI"])]
if mi:print("  缺数据集:",mi)
print("  mean ARI=%.3f  NMI=%.3f  ACC=%.3f"%(
 np.nanmean([per[d]["ARI"] for d in LINEUP]),np.nanmean([per[d]["NMI"] for d in LINEUP]),np.nanmean([per[d]["CA"] for d in LINEUP])))
for d in LINEUP:print("    %-24s ARI=%.3f NMI=%.3f ACC=%.3f"%(d,per[d]["ARI"],per[d]["NMI"],per[d]["CA"]))
print("="*62);print("CHECK2  消融 15-set  Δ = full - 变体");print("="*62)
ABL=f"{PROJ}/results/ablation"
ds=sorted(d for d in (os.path.basename(x.rstrip('/')) for x in glob.glob(f"{ABL}/full/*/")) if d not in EX)
print("  消融数据集数(应为15):",len(ds))
full={x:rd(f"{ABL}/full/{x}") for x in ds}
print("  组件        ARI       NMI       ACC      3指标平均")
for var in ["wo_cl","wo_graph","wo_edge","wo_recon"]:
    vm={x:rd(f"{ABL}/{var}/{x}") for x in ds}
    o={mk:np.nanmean([full[x][mk]-vm[x][mk] for x in ds]) for mk in ["ARI","NMI","CA"]}
    print("  %-9s %+.4f  %+.4f  %+.4f   %+.4f"%(var,o["ARI"],o["NMI"],o["CA"],np.nanmean([o["ARI"],o["NMI"],o["CA"]])))
print("  → 这组 15-set 数字填回 .tex 消融段")
