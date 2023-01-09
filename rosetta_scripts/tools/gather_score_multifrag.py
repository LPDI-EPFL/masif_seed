import os, sys
import pandas as pd
import re


cmd = 'ls ./results | wc -l'
nseed= int(os.popen(cmd).readlines()[0])
df=pd.DataFrame([])
for seedid in range(1,nseed+1,1):
    for array in range(1,363):
        path='./results/seed_'+str(seedid)+'/'+str(array)+'/score_S'+str(seedid)+'.sc'
        if (os.path.exists(path)):
            f = open(path,'r')
            outname='tmp_'+str(seedid)+'_'+str(array)+'.sc'
            tmpf = open(outname,'w')
            for line in f:
                metrics=re.split(' +',line.strip())
                new_metrics=metrics.copy()
                if metrics[0]=='SCORE:':
                    for element in metrics:
                        if (element=='0') or (',' in element):
                            new_metrics.remove(element)
                        if (element in ['graft_in_motif_ranges','graft_in_scaffold_ranges', 'graft_out_scaffold_ranges', 'graft_scaffold_size_change']):
                            new_metrics.remove(element)
                for i in range(0,len(new_metrics),1):
                    tmpf.write(new_metrics[i])
                    if i==(len(new_metrics)-1):
                        tmpf.write('\n')
                    else:
                        tmpf.write(';')
            f.close()
            tmpf.close()
            df1 = pd.read_csv(outname, sep=';',skiprows=1).dropna()

            df1 = df1[df1.total_score != 'total_score'].drop(columns=['SCORE:'])
            df1 = df1.applymap(lambda x : pd.to_numeric(x,errors='ignore',downcast='float'))
            df1['seed']=str(seedid)
            df1['array']=str(array)
            df1.shape
            df = pd.concat([df,df1])

            os.remove(outname)
#df = df.drop(columns=['dslf_fa13','fa_atr','fa_dun','fa_elec','fa_intra_rep','fa_intra_sol_xover4','fa_rep','fa_sol','graft_max_motif_fragment_RMSD','hbond_bb_sc','hbond_lr_bb','hbond_sc','hbond_sr_bb','lk_ball_wtd','omega','p_aa_pp','pro_close','rama_prepro','ref','yhh_planarity','graft_full_bb_mode','graft_in_motif_ranges','graft_in_scaffold_ranges','graft_out_scaffold_ranges','graft_scaffold_size_change'])
df=df[['total_score','buried_unsat_hbonds','sasa','sc','hbonds','ddG','graft_RMSD','array','seed','description']]
#print(df)
for column in list(df):
    if df[column].dtypes=='float64':
        df[column] = df[column].map('{:.3f}'.format)
df.reset_index(inplace=True)
dico_ori={}
for idx in df.index:
    desc=str(df.at[idx, 'description'])
    ori='NA'
    if(desc[0:3]=='AF-'):
        ori='af'
    elif((desc[0:4]=='HHH_') or (desc[0:5]=='EHEE_') or (desc[0:6]=='EEHEE_') or (desc[0:6]=='EEHEE_') or (desc[0:5]=='HEEH_') or (desc[0:4]=='HHH_') or (desc[0:3]=='HH_') or (desc[0:5]=='HEEE_') or (desc[0:4]=='EHE_') or (desc[0:5]=='EEHE_') or (desc[0:4]=='EEH_') or (desc[0:5]=='EEEH_') or (desc[0:7]=='EEE_EEE') or (desc[0]=='g')):
        ori='miniprotein'
    else:
        ori='pdb'
    dico_ori[desc]=ori
df['origin'] = df['description'].map(dico_ori)

df.to_csv('all_score.sc',sep='\t',index=False)
