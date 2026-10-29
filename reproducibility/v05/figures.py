"""Publication figures for primary provenance, deduction and external transfer."""
from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'reviews/v05';FIG=ROOT/'paper/figures'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.dpi':180})
COL=['#196B75','#D27A36','#72618C','#677B40']
def save(fig,name):
 FIG.mkdir(exist_ok=True);fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight',metadata={'CreationDate':None,'ModDate':None});fig.savefig(OUT/(name+'.png'),bbox_inches='tight');plt.close(fig)
def label(ax,text):ax.set_title(text,loc='left',fontweight='bold',pad=10)

def primary():
 sites=pd.read_csv(OUT/'primary_sites.tsv',sep='\t')
 fig=plt.figure(figsize=(9,6.4),layout='constrained');gs=fig.add_gridspec(2,2,height_ratios=[1.25,1])
 ax=fig.add_subplot(gs[0,:]);pts=ax.scatter(sites.longitude,sites.latitude,c=sites.campaigns,s=25+sites.samples*1.2,cmap='viridis',vmin=1,vmax=5,edgecolor='white',linewidth=.6)
 aux=sites[sites.site.gt(60)];ax.scatter(aux.longitude,aux.latitude,s=75,facecolor='none',edgecolor='#9B452F',linewidth=1.3)
 ax.set(xlabel='Longitude (°E)',ylabel='Latitude (°N)',ylim=(18.95,21.15));label(ax,'A  Primary sampling footprint in the Rub’ al-Khali')
 cb=fig.colorbar(pts,ax=ax,shrink=.65,pad=.015,ticks=[1,2,3,4,5]);cb.set_label('Campaigns sampled')
 ax.text(.01,.03,'551 samples · 61 sites · five campaigns\nPoint area: sample count; outlined points: auxiliary sites',transform=ax.transAxes,fontsize=8)
 ax=fig.add_subplot(gs[1,0]);panel=pd.read_csv(ROOT/'data/v03/sample_evidence.tsv',sep='\t');counts=pd.crosstab(panel.trip_num,panel.compartment)
 bottom=np.zeros(5)
 for i,c in enumerate(['Surface','Deep','PlantRoot']):ax.bar(np.arange(5),counts[c],bottom=bottom,label={'PlantRoot':'Root-adjacent','Deep':'Subsurface','Surface':'Surface'}[c],color=COL[i]);bottom+=counts[c].to_numpy()
 ax.set(xticks=np.arange(5),xticklabels=['Mar ’23','Jul ’23','Feb ’24','Aug ’24','Oct ’25'],ylabel='Metagenomes');ax.tick_params(axis='x',rotation=25);ax.legend(frameon=False,fontsize=8);label(ax,'B  Repeated biological observations')
 ax=fig.add_subplot(gs[1,1]);ax.axis('off');label(ax,'C  From observations to intervention effects')
 boxes=[(.5,.9,'Field specimens → sequences'),(.5,.65,'Genomes → deduced pathway support'),(.5,.4,'Admitted reactions + physical mechanisms'),(.5,.15,'Specified intervention → conditional effect')]
 for i,(x,y,t) in enumerate(boxes):
  ax.add_patch(FancyBboxPatch((.01,y-.075),.98,.15,boxstyle='round,pad=.015',facecolor='#EFF4F3' if i<2 else '#F8F1E9',edgecolor='#BECAC8',transform=ax.transAxes));ax.text(x,y,t,ha='center',va='center',transform=ax.transAxes,fontsize=8.5)
  if i<3:ax.annotate('',xy=(.5,y-.17),xytext=(.5,y-.085),xycoords='axes fraction',arrowprops={'arrowstyle':'->','color':'#657271'})
 save(fig,'v05_primary_twin')

def deduction():
 d=pd.read_csv(OUT/'capacity_overlap_bounds.tsv',sep='\t').dropna().sort_values('upper');r=json.loads((OUT/'deduction_causal_trace.json').read_text())
 fig,axes=plt.subplots(1,2,figsize=(9,3.6),layout='constrained')
 ax=axes[0];x=np.arange(len(d));ax.fill_between(x,d.lower*100,d.upper*100,color='#DDD9CC',label='Marginal-only admissible range');ax.scatter(x,d.genome_resolved*100,s=3,c=COL[0],label='Genome-resolved intersection')
 ax.set(xlabel='Samples ordered by overlap upper bound',ylabel='H₂ + CBB carrier share (%)');ax.legend(frameon=False,fontsize=7.5,loc='upper left');label(ax,'A  Deduction resolves an ambiguous conjunction')
 ax=axes[1];keys=['ontology_supported','same_information_flat','unsupported_CO_admission'];vals=[r[k]['fluxes']['carbohydrate_c']*1e9 for k in keys]
 ax.barh([2,1,0],vals,color=[COL[0],COL[0],COL[1]],height=.55);ax.set(yticks=[2,1,0],yticklabels=['Ontology + admission','Same-information flat','Unsupported CO admission'],xlabel='Conditional C ceiling (nmol m⁻² s⁻¹)',xlim=(0,.59));label(ax,'B  Reaction admission changes the computed bound')
 for y,v in zip([2,1,0],vals):ax.text(v+.008,y,f'{v:.3f}',va='center',fontsize=8)
 save(fig,'v05_deduction')

def transfer():
 report=json.loads((OUT/'external_transfer.json').read_text());nam=pd.read_csv(OUT/'namib_observations.tsv',sep='\t');lab=pd.read_csv(OUT/'sonoran_lab_observations.tsv',sep='\t')
 fig,axes=plt.subplots(2,2,figsize=(9,7),layout='constrained')
 ax=axes[0,0]
 for i,site in enumerate(['C14-2','C14-4']):
  means=[]
  for j,state in enumerate(['dry','wet']):
   y=nam[(nam.site==site)&nam.state.eq(state)].nmol_h_g.to_numpy();means.append(y.mean());ax.scatter(np.repeat(j+i*.05,len(y)),y,c=COL[i],s=15,alpha=.65)
  ax.plot(np.array([0,1])+i*.05,means,color=COL[i],label=site,marker='o')
 ax.set(xticks=[0,1],xticklabels=['Native dry','Water added'],yscale='log',ylabel='H₂ uptake (nmol h⁻¹ g⁻¹)');ax.legend(frameon=False);label(ax,'A  Namib hydration response')
 ax=axes[0,1]
 for i,(name,g) in enumerate(lab.groupby('Treatment')):
  means=g.groupby('gwc').uptake_nmol_m2_s.mean();ax.plot(means.index,means.values,'o-',label=name,color=COL[i],lw=1.4,ms=4)
 ax.set(xlabel='Gravimetric water content (g g⁻¹)',ylabel='H₂ uptake (nmol m⁻² s⁻¹)',xticks=[.02,.10,.20]);ax.axhline(0,color='#aaa',lw=.6);ax.legend(frameon=False,fontsize=7);label(ax,'B  Sonoran controlled moisture response')
 ax=axes[1,0];data={r['model']:r['rmse'] for r in report['sonoran']['field_metrics']}
 names=['frozen_reference','calibrated_kinetics','calibrated_diffusion_only','linear_moisture_temperature','training_mean','calibrated_retention'];labels=['Frozen reference','Kinetics calibration','Diffusion-only calibration','Linear moisture + temperature','Training mean','Soil-response calibration']
 ax.barh(np.arange(6),[data[n] for n in names],color=['#B5BDBB',COL[1],COL[1],'#8A9092','#8A9092',COL[0]],height=.65);ax.set(yticks=np.arange(6),yticklabels=labels,xlabel='Held-campaign RMSE (nmol m⁻² s⁻¹)',xlim=(0,2.1));ax.invert_yaxis();ax.tick_params(axis='y',labelsize=8);label(ax,'C  Sonoran field transfer requires calibration')
 for y,n in enumerate(names):ax.text(data[n]+.025,y,f'{data[n]:.2f}',va='center',fontsize=8)
 ax=axes[1,1];d=pd.DataFrame(report['temperate_leave_soil_out']);soils=['Forest','Meadow','PB_sand'];models=['diffusion_only','training_mean','activation_diffusion']
 for j,m in enumerate(models):ax.bar(np.arange(3)+(j-1)*.23,[d[(d.soil==s)&d.model.eq(m)].rmse.iloc[0] for s in soils],width=.21,label={'diffusion_only':'Diffusion','training_mean':'Training mean','activation_diffusion':'Activation + diffusion'}[m],color=['#AEB7B5','#8A9092',COL[0]][j])
 ax.set(xticks=np.arange(3),xticklabels=['Forest','Meadow','Sand'],ylabel='Held-soil bottle-rate RMSE (min⁻¹)');ax.legend(frameon=False,fontsize=7);label(ax,'D  Cross-soil transfer varies by environment')
 save(fig,'v05_external_transfer')
if __name__=='__main__':primary();deduction();transfer()
