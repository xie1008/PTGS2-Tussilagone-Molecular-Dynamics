import os,pickle,csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.vq import kmeans2
ROOT=os.path.abspath(os.path.join(os.path.dirname(__file__),'../..')); OUT=os.path.join(ROOT,'outputs/PTGS2_TS_Selected_3Rep_Final_QC3')
with open('work/rep4_all_data.pkl','rb') as f: all_data,atoms,ref,ca,prot,lig,pock=pickle.load(f)
ORIG=[2,3,4]; DISP={2:1,3:2,4:3}; COL={1:'#1F77B4',2:'#FF7F0E',3:'#2CA02C'}
data={DISP[r]:all_data[r] for r in ORIG}
def rms(a,b): return float(np.sqrt(np.mean(np.sum((a-b)**2,axis=1))))
def wr(name,head,rows):
 with open(os.path.join(OUT,'CSV',name),'w',newline='') as f: w=csv.writer(f);w.writerow(head);w.writerows(rows)
with open(os.path.join(OUT,'provenance_mapping.csv'),'w',newline='') as f:
 w=csv.writer(f);w.writerow(['display_replica','original_replica','seed','platform_precision','source_directory','selection_status']);
 w.writerow([1,2,20261107,'Apple M4/OpenCL single','work/ptgs2_rep2','selected for this reanalysis']);w.writerow([2,3,20261108,'Tesla T4/OpenCL mixed','work/ptgs2_rep3','selected for this reanalysis']);w.writerow([3,4,20261109,'Tesla T4/OpenCL mixed','work/ptgs2_rep4','selected for this reanalysis'])
summary=[];lateX=[];meta=[]
for r,d in data.items():
 idx=int(np.argmin(abs(d['t']-50))); rr=np.array([rms(x,d['ligcoords'][idx]) for x in d['ligcoords']]);
 wr(f'replica{r}_ligand_rmsd.csv',['time_ns','rmsd_A','reference_time_ns','reference_frame_index','original_replica'],[[t,v,50,idx,ORIG[r-1]] for t,v in zip(d['t'],rr)])
 wr(f'replica{r}_contact_occupancy.csv',['chain','residue_number','residue','occupancy_percent','original_replica'],[[*x,ORIG[r-1]] for x in d['occ']])
 late=rr[d['t']>=50];w1=rr[(d['t']>=50)&(d['t']<75)];w2=rr[d['t']>=75]
 for i,t in enumerate(d['t']):
  if t>=50:lateX.append(d['ligcoords'][i].reshape(-1));meta.append((r,float(t)))
 summary.append([r,ORIG[r-1],d['protein_r'][d['t']>=50].mean(),d['rg'][d['t']>=50].mean(),late.mean(),late.std(ddof=1),w1.mean(),w2.mean(),d['contacts'][d['t']>=50].mean(),d['mind'][d['t']>=50].mean()])
X=np.asarray(lateX);cent,lab=kmeans2(X,3,minit='++',iter=30,seed=7);order=np.argsort([sum(lab==k) for k in range(3)])[::-1];mp={old:new for new,old in enumerate(order)};lab=np.array([mp[x] for x in lab]);cent=cent[order]
wr('joint_cluster_assignments.csv',['display_replica','time_ns','cluster','original_replica'],[[r,t,int(c+1),ORIG[r-1]] for (r,t),c in zip(meta,lab)])
crow=[]
for r in range(1,4):
 q=[c for (rr,t),c in zip(meta,lab) if rr==r]
 for c in range(3):crow.append([r,c+1,sum(x==c for x in q),100*sum(x==c for x in q)/len(q)])
wr('joint_cluster_occupancy.csv',['display_replica','cluster','frames','percent','original_replica'],[[r,c,n,p,ORIG[r-1]] for r,c,n,p in crow])
repmed=[]
for r in range(1,4):
 ix=np.array([q for q,(rr,t) in enumerate(meta) if rr==r]);center=X[ix].mean(0);q=ix[np.argmin(np.sum((X[ix]-center)**2,axis=1))];repmed.append(X[q].reshape(28,3))
K=np.array([[rms(a,b) for b in repmed] for a in repmed]);wr('cross_replica_medoid_rmsd_3x3.csv',['display_row','display_col','rmsd_A'],[[i+1,j+1,K[i,j]] for i in range(3) for j in range(3)])
occ={(c,n,res,r):pct for r,d in data.items() for c,n,res,pct in d['occ']};keys=sorted(set((c,n,res) for c,n,res,r in occ if sum(occ.get((c,n,res,rr),0)>=50 for rr in range(1,4))>=2),key=lambda z:(z[0],z[1]));M=np.array([[occ.get((c,n,res,r),0) for r in range(1,4)] for c,n,res in keys]);wr('key_contact_occupancy.csv',['chain','residue_number','residue','Replica1_percent','Replica2_percent','Replica3_percent'],[[c,n,res,*M[i]] for i,(c,n,res) in enumerate(keys)])
top={r:', '.join(f'{c}{n}{res} ({pct:.1f}%)' for c,n,res,pct in data[r]['occ'][:6]) for r in data};wr('summary_table.csv',['display_replica','original_replica','protein_RMSD_last50_A','Rg_last50_A','ligand_RMSD_last50_mean_A','ligand_RMSD_last50_SD_A','ligand_RMSD_50_75_mean_A','ligand_RMSD_75_100_mean_A','contacts_last50_mean','min_distance_last50_A','major_contacts'],[[*s,top[s[0]]] for s in summary])
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.linewidth':.8,'pdf.fonttype':42,'ps.fonttype':42})
def export(fig,folder,name):
 os.makedirs(folder,exist_ok=True)
 for ext,kw in [('pdf',{}),('svg',{}),('png',{'dpi':300}),('tiff',{'dpi':300})]:fig.savefig(os.path.join(folder,name+'.'+ext),bbox_inches='tight',pad_inches=.05,**kw)
 plt.close(fig)
def one(ylabel,xlabel,series,legend=True):
 f,a=plt.subplots(figsize=(3.25,2.55),constrained_layout=True)
 for r,x in series:a.plot(x[0],x[1],color=COL[r],lw=.85,label=f'Replica {r}')
 a.set(xlabel=xlabel,ylabel=ylabel);a.grid(axis='y',color='.9',lw=.5)
 if legend:a.legend(frameon=False,loc='best')
 return f
fig,ax=plt.subplots(1,3,figsize=(10.2,2.9),constrained_layout=True)
for r,d in data.items():
 idx=int(np.argmin(abs(d['t']-50)));rr=np.array([rms(x,d['ligcoords'][idx]) for x in d['ligcoords']]);ax[0].plot(d['t'],d['protein_r'],color=COL[r],label=f'Replica {r}',lw=.85);ax[1].plot(d['t'],rr,color=COL[r],label=f'Replica {r}',lw=.85);ax[2].plot(d['t'],d['contacts'],color=COL[r],label=f'Replica {r}',lw=.85)
for a,y in zip(ax,['Protein Cα RMSD (Å)','Pocket-aligned ligand RMSD (Å)','Heavy-atom contacts (≤4.5 Å)']):a.set(xlabel='Time (ns)',ylabel=y,xlim=(0,100));a.grid(axis='y',color='.9',lw=.5)
ax[1].text(.02,.04,'fit: Chain A A43–A89\nreference: own 50 ns frame',transform=ax[1].transAxes,fontsize=7);ax[2].legend(frameon=False,ncol=3,loc='upper center',bbox_to_anchor=(.5,-.25));export(fig,os.path.join(OUT,'Combined_Figures'),'Figure_8_combined')
series1=[(r,(d['t'],d['protein_r'])) for r,d in data.items()]
series2=[]
series3=[(r,(d['t'],d['contacts'])) for r,d in data.items()]
for r,d in data.items():
 idx=int(np.argmin(abs(d['t']-50)));series2.append((r,(d['t'],[rms(x,d['ligcoords'][idx]) for x in d['ligcoords']])))
for k,(y,series) in enumerate([('Protein Cα RMSD (Å)',series1),('Pocket-aligned ligand RMSD (Å)',series2),('Heavy-atom contacts (≤4.5 Å)',series3)],1):
 export(one(y,'Time (ns)',series),os.path.join(OUT,'Individual_Panels',f'Figure_8{k}'),'Figure_8'+chr(64+k))

# Supplementary six-panel figure
fig,ax=plt.subplots(2,3,figsize=(10.2,6.1),constrained_layout=True);ax=ax.ravel()
for r,d in data.items():
 idx=int(np.argmin(abs(d['t']-50)));rr=np.array([rms(x,d['ligcoords'][idx]) for x in d['ligcoords']]);ax[0].plot(np.arange(len(d['rmsf'])),d['rmsf'],color=COL[r],label=f'Replica {r}',lw=.8);ax[1].plot(d['t'],d['rg'],color=COL[r],label=f'Replica {r}',lw=.8)
for a,y,x in [(ax[0],'Protein Cα RMSF (Å)','Model residue index'),(ax[1],'Radius of gyration (Å)','Time (ns)')]:a.set(xlabel=x,ylabel=y);a.legend(frameon=False,fontsize=7)
im=ax[2].imshow(M,aspect='auto',cmap='viridis',vmin=0,vmax=100);ax[2].set(xticks=[0,1,2],xticklabels=['Replica 1','Replica 2','Replica 3'],xlabel='Replica',ylabel='Model residue');ax[2].set_yticks(range(len(keys)),[f'{c}:{n} {res}' for c,n,res in keys],fontsize=7);fig.colorbar(im,ax=ax[2],label='Occupancy (%)')
bottom=np.zeros(3)
for c in range(3):
 p=np.array([next(x[3] for x in crow if x[0]==r and x[1]==c+1) for r in range(1,4)]);ax[3].bar(['Replica 1','Replica 2','Replica 3'],p,bottom=bottom,label=f'C{c+1}');bottom+=p
ax[3].set(xlabel='Replica',ylabel='Cluster occupancy (%)');ax[3].legend(frameon=False,ncol=3,fontsize=7,loc='lower center',bbox_to_anchor=(.5,1.02),borderaxespad=0)
for r,d in data.items():
 idx=int(np.argmin(abs(d['t']-50)));rr=np.array([rms(x,d['ligcoords'][idx]) for x in d['ligcoords']]);ax[4].hist(rr[d['t']>=50],bins=np.arange(0,12.01,.4),density=True,histtype='step',lw=1.2,color=COL[r],label=f'Replica {r}')
ax[4].set(xlabel='Ligand RMSD (Å)',ylabel='Probability density');ax[4].legend(frameon=False,fontsize=7)
hm=ax[5].imshow(K,cmap='magma',vmin=0,vmax=max(1,K.max()));ax[5].set(xticks=[0,1,2],xticklabels=['Replica 1','Replica 2','Replica 3'],yticks=[0,1,2],yticklabels=['Replica 1','Replica 2','Replica 3'],xlabel='Representative pose',ylabel='Representative pose')
for i in range(3):
 for j in range(3):
  rgba=plt.get_cmap('magma')(K[i,j]/max(1,K.max()));lum=0.2126*rgba[0]+0.7152*rgba[1]+0.0722*rgba[2];ax[5].text(j,i,f'{K[i,j]:.2f}',ha='center',va='center',color='black' if lum>0.55 else 'white')
fig.colorbar(hm,ax=ax[5],label='RMSD (Å)');export(fig,os.path.join(OUT,'Combined_Figures'),'Figure_S8_combined')
for letter,ylabel in zip('ABCDEF',['Protein Cα RMSF (Å)','Radius of gyration (Å)','Contact occupancy (%)','Cluster occupancy (%)','Probability density','RMSD (Å)']):
 f,a=plt.subplots(figsize=(3.25,2.55),constrained_layout=True)
 if letter=='A':
  for r,d in data.items():a.plot(np.arange(len(d['rmsf'])),d['rmsf'],color=COL[r],label=f'Replica {r}',lw=.8)
  a.set(xlabel='Model residue index',ylabel=ylabel);a.legend(frameon=False,fontsize=7)
 elif letter=='B':
  for r,d in data.items():a.plot(d['t'],d['rg'],color=COL[r],label=f'Replica {r}',lw=.8)
  a.set(xlabel='Time (ns)',ylabel=ylabel);a.legend(frameon=False,fontsize=7)
 elif letter=='C':
  a.imshow(M,aspect='auto',cmap='viridis',vmin=0,vmax=100);a.set(xticks=[0,1,2],xticklabels=['R1','R2','R3'],xlabel='Replica',ylabel='Model residue');a.set_yticks(range(len(keys)),[f'{c}:{n} {res}' for c,n,res in keys],fontsize=7)
 elif letter=='D':
  b=np.zeros(3)
  for c in range(3):
   p=np.array([next(x[3] for x in crow if x[0]==r and x[1]==c+1) for r in range(1,4)]);a.bar(['R1','R2','R3'],p,bottom=b,label=f'C{c+1}');b+=p
  a.set(xlabel='Replica',ylabel=ylabel);a.legend(frameon=False,ncol=3,fontsize=7)
 elif letter=='E':
  for r,d in data.items():
   idx=int(np.argmin(abs(d['t']-50)));rr=np.array([rms(x,d['ligcoords'][idx]) for x in d['ligcoords']]);a.hist(rr[d['t']>=50],bins=np.arange(0,12.01,.4),density=True,histtype='step',lw=1.2,color=COL[r],label=f'Replica {r}')
  a.set(xlabel='Ligand RMSD (Å)',ylabel=ylabel);a.legend(frameon=False,fontsize=7)
 else:
  a.imshow(K,cmap='magma',vmin=0,vmax=max(1,K.max()));a.set(xticks=[0,1,2],xticklabels=['R1','R2','R3'],yticks=[0,1,2],yticklabels=['R1','R2','R3'],xlabel='Representative pose',ylabel='Representative pose')
  for i in range(3):
   for j in range(3):a.text(j,i,f'{K[i,j]:.2f}',ha='center',va='center',color='white' if K[i,j]>K.max()/2 else 'black')
 export(f,os.path.join(OUT,'Individual_Panels',f'Figure_S8{letter}'),f'Figure_S8{letter}')
