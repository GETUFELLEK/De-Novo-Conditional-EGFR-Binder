"""Rigid independent-site proton-linkage approximation; not affinity or a validated switch."""
import argparse
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser,PDBIO
from geometry import dump
try:
 from propka.run import single
except ImportError:
 single=None

def linkage(path,out,binder_chains,target_chains):
 out.mkdir(parents=True,exist_ok=True);m=PDBParser(QUIET=True).get_structure('p',path)[0];states={}
 binder=''.join(binder_chains);target=''.join(target_chains)
 assert not set(binder)&set(target)
 assert set(c.id for c in m)==set(binder+target), 'Assign every present chain explicitly'
 for label,chains in [('complex',target+binder),('binder',binder),('target',target)]:
  sub=m.copy()
  for c in list(sub):
   if c.id not in chains:sub.detach_child(c.id)
  fp=out/(label+'.pdb');io=PDBIO();io.set_structure(sub);io.save(str(fp));mol=single(str(fp),write_pka=False)
  states[label]={f'{g.atom.chain_id}:{g.atom.res_num}:{g.residue_type}':dict(pka=float(g.pka_value),buried=float(g.buried),coupled=bool(g.non_covalently_coupled_groups or g.covalently_coupled_groups)) for g in mol.conformations['AVR'].groups if g.titratable and not g.exclude_cys_from_results}
 bound=states['complex'];free=dict(states['binder'],**states['target']);assert set(bound)==set(free)
 bs=np.array([bound[k]['pka'] for k in bound]);fs=np.array([free[k]['pka'] for k in bound]);rt=.00198720425864083*298.15
 def energy(pH):return float(-rt*np.sum(np.logaddexp(0,np.log(10)*(bs-pH))-np.logaddexp(0,np.log(10)*(fs-pH))))
 x=dict(deltaG_6_minus7p4_kcal_mol=energy(6)-energy(7.4),bound=bound,free=free,
  interpretation='Rigid independent-site proton-linkage approximation; not total binding free energy or measured selectivity')
 dump(out/'analysis.json',x);return x

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--binder-chains',nargs='+',required=True);p.add_argument('--target-chains',nargs='+',required=True);a=p.parse_args()
 if single is None:p.error('Install PROPKA in the analysis environment')
 linkage(a.input,a.output,a.binder_chains,a.target_chains)
