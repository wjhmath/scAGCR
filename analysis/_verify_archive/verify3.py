import numpy as np, glob, os, json
from scipy.optimize import linear_sum_assignment
try:
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
except Exception:
    adjusted_rand_score=None
PROJ="/home/liyang/BioJiaheWang/scGTAC"

print("="*62);print("CHECK3  Muraro 案例 ARI + 逐类型 recall (Fig3 源)");print("="*62)
cand=glob.glob(f"{PROJ}/paper_figures/Figure3_umap/emb_muraro_scGTAC.npz")+glob.glob(f"{PROJ}/paper_figures*/*muraro*scGTAC*.npz")
if not cand:
    print("  找不到 npz,改路径")
else:
    z=np.load(cand[0],allow_pickle=True)
    y=list(map(str,np.asarray(z['y']).ravel())); pred=list(map(str,np.asarray(z['pred']).ravel()))
    if adjusted_rand_score:
        print("  ARI(重算)=%.3f   NMI=%.3f   (论文案例 ARI=0.913)"%(
            adjusted_rand_score(y,pred), normalized_mutual_info_score(y,pred)))
    types=sorted(set(y)); clus=sorted(set(pred))
    C=np.zeros((len(types),len(clus)))
    for a,b in zip(y,pred): C[types.index(a),clus.index(b)]+=1
    r,cc=linear_sum_assignment(-C)
    rec={t:0.0 for t in types}
    for i,j in zip(r,cc):
        if C[i].sum()>0: rec[types[i]]=C[i,j]/C[i].sum()
    print("  逐类型 recall (Hungarian 1-1, 对照 Fig3b):")
    for t in types:
        print("    %-30s n=%-5d recall=%.3f"%(t,int(C[types.index(t)].sum()),rec[t]))
    print("  论文写的: alpha 0.97, delta 0.97, PP 1.00, mesen 1.00, beta 0.95, endoth 0.95, ductal 0.87, acinar 0.59")

print("="*62);print("CHECK5  baseline JSON 格式 (看键名)");print("="*62)
js=glob.glob(f"{PROJ}/results/baselines/seurat/muraro_pancreas_seed1_metrics.json")+\
   glob.glob(f"{PROJ}/results/baselines/*/muraro_pancreas_seed1_metrics.json")
if js:
    print("  示例文件:",js[0]); print("  内容:",json.load(open(js[0])))
else:
    print("  没找到,贴: find results/baselines* -name '*metrics.json' | head")
print("  baseline 方法目录:")
for base in ["results/baselines","results/baselines_chuli"]:
    p=f"{PROJ}/{base}"
    if os.path.isdir(p): print("   ",base,"->",sorted(os.listdir(p)))
