"""Fixed-backbone TNF/Fv side-chain packing; energies are model scores, not affinities."""
import argparse,json
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser,PDBIO,Select
from Bio.SeqUtils import seq3
import pyrosetta as py
from pyrosetta import rosetta as ro
from geometry import REGIONS,dump
class Heavy(Select):
 def accept_atom(self,a):return a.element not in ['H','D']
def main(a):
 py.init(f'-mute all -load_PDB_components false -ex1 -ex2aro -use_input_sc -constant_seed -jran {a.seed}')
 p=py.pose_from_pdb(str(a.input));orig=p.clone();info=p.pdb_info();sf=py.create_score_function('ref2015');region=({('D',p.pdb_info().number(i)) for i in range(1,p.size()+1) if p.pdb_info().chain(i)=='D'} if a.molecule_class=='compact' else {(k[0],n) for k,v in REGIONS.items() for n in v})
 if a.sequences:
  seq=json.loads(a.sequences.read_text())
  for c in ('D' if a.molecule_class=='compact' else 'HL'):
   for n,aa in enumerate(seq[c],1):
    i=info.pdb2pose(c,n);assert i
    if p.residue(i).name1()!=aa:
     assert (c,n) in region;mut=ro.protocols.simple_moves.MutateResidue(i,seq3(aa).upper());mut.set_preserve_atom_coords(True);mut.apply(p)
 def metrics(p):
  q=p.clone();ia=ro.protocols.analysis.InterfaceAnalyzerMover('D_ABC' if a.molecule_class=='compact' else 'HL_ABC');ia.set_scorefunction(sf);ia.set_pack_input(False);ia.set_pack_separated(False);ia.set_compute_interface_sc(True);ia.apply(q);return {str(k):float(v) for k,v in q.scores.items() if isinstance(v,(float,int))}
 before=metrics(p);centers=np.array([list(p.residue(info.pdb2pose(c,n)).xyz('CA')) for c,n in region]);mov=[]
 for i in range(1,p.size()+1):
  c,n=info.chain(i),info.number(i)
  if (c,n) in region:mov.append(i)
 if a.repack_target_interface:
  binder_atoms=np.array([list(p.residue(i).xyz(j)) for i in mov for j in range(1,p.residue(i).nheavyatoms()+1)])
  from scipy.spatial import cKDTree
  tree=cKDTree(binder_atoms)
  for i in range(1,p.size()+1):
   if info.chain(i) in 'ABC' and min(tree.query(np.array([list(p.residue(i).xyz(j)) for j in range(1,p.residue(i).nheavyatoms()+1)]))[0])<=6:mov.append(i)
 task=py.standard_packer_task(p);task.restrict_to_repacking();task.or_include_current(True);mm=ro.core.kinematics.MoveMap();mm.set_bb(False);mm.set_chi(False);mm.set_jump(False)
 for i in range(1,p.size()+1):
  if i not in mov:task.nonconst_residue_task(i).prevent_repacking()
  else:mm.set_chi(i,True)
 ro.protocols.minimization_packing.PackRotamersMover(sf,task).apply(p);mini=ro.protocols.minimization_packing.MinMover(mm,sf,'lbfgs_armijo_nonmonotone',.001,True);mini.max_iter(150);mini.apply(p)
 shift=max(float(np.linalg.norm(np.array(list(p.residue(i).xyz(at)))-np.array(list(orig.residue(i).xyz(at))))) for i in range(1,p.size()+1) for at in ['N','CA','C','O']);assert shift<1e-4,shift
 fixed_backbone_metrics=metrics(p);movement={};bounded=True
 if a.sequences:
  cart=py.create_score_function('ref2015_cart');cart.set_weight(ro.core.scoring.coordinate_constraint,1)
  for i in range(1,p.size()+1):
   key=info.chain(i),info.number(i);res=p.residue(i)
   if key in region:mm.set_bb(i,True)
   for at in ['N','CA','C','O']:
    p.add_constraint(ro.core.scoring.constraints.CoordinateConstraint(ro.core.id.AtomID(res.atom_index(at),i),ro.core.id.AtomID(2,1),res.xyz(at),ro.core.scoring.func.HarmonicFunc(0,.35 if key in region else .05)))
  local=ro.protocols.minimization_packing.MinMover(mm,cart,'lbfgs_armijo_nonmonotone',.001,True);local.cartesian(True);local.max_iter(200);local.apply(p);p.remove_constraints()
  for group in ['CDR','fixed']:
   values=[float(np.linalg.norm(np.array(list(p.residue(i).xyz(at)))-np.array(list(orig.residue(i).xyz(at))))) for i in range(1,p.size()+1) if ((info.chain(i),info.number(i)) in region)==(group=='CDR') for at in ['N','CA','C','O']];movement[group]=dict(RMSD_A=float(np.sqrt(np.mean(np.square(values)))),max_A=max(values))
  bounded=movement['CDR']['RMSD_A']<=.5 and movement['CDR']['max_A']<=1.5 and movement['fixed']['max_A']<=.01
 after=metrics(p);a.output.mkdir(parents=True,exist_ok=True);p.dump_pdb(str(a.output/'rosetta.pdb'));m=PDBParser(QUIET=True).get_structure('p',a.output/'rosetta.pdb')[0];io=PDBIO();io.set_structure(m);io.save(str(a.output/'protein.pdb'),Heavy())
 # Predefined interface quality gate is applied together with post-packing physical validation.
 gate=bounded and after['sc_value']>=.5 and after['dSASA_int']>=500 and after['hbonds_int']>=2 and after['dG_separated']<0
 dump(a.output/'packing.json',dict(protocol_version=4 if a.repack_target_interface else 3,sidechain_scope='Binder and target interface chi only; all backbones fixed' if a.repack_target_interface else 'CDR chi only, then sequence-specific restrained CDR Cartesian minimization; target/framework fixed',target_interface_repack=a.repack_target_interface,movable_residues=[(info.chain(i),info.number(i)) for i in mov],fixed_backbone_metrics=fixed_backbone_metrics,movement=movement,bounded_movement=bounded,input=str(a.input),output=str(a.output/'protein.pdb'),before=before,after=after,backbone_max_shift_A=shift,packing_pass=gate,seed=a.seed,interpretation='Rosetta interface energy supports local packing; not experimental binding energy or pH selectivity',thresholds=dict(sc_min=.5,buried_area_min_A2=500,interface_hbond_min=2,interface_REU_max=0)))
 print('PACK',gate,{k:after[k] for k in ['sc_value','dSASA_int','hbonds_int','dG_separated']},flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--sequences',type=Path);p.add_argument('--molecule-class',choices=['Fab','compact'],default='Fab');p.add_argument('--seed',type=int,default=66000);p.add_argument('--repack-target-interface',action='store_true');a=p.parse_args();assert not(a.sequences and a.repack_target_interface),'Separate sequence changes from target-rotamer rescue';main(a)
