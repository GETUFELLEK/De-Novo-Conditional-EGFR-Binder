"""Explicit multichain TNF/Fv geometry; distances in angstroms."""
import json,os
PILOT=os.environ.get("TNF_PILOT","pilot01")
H3_LENGTH=int(os.environ.get("TNF_H3_LENGTH","17"))
TARGET_START=int(os.environ.get("TNF_TARGET_START","6"))
TARGET_FILE="human_tnf_canonical_core.pdb" if TARGET_START==6 else "human_tnf_canonical_full.pdb"
import numpy as np
from scipy.spatial import cKDTree
from Bio.PDB.vectors import calc_angle,calc_dihedral
EXPECTED={aa:set(('N CA C O '+sc).split()) for aa,sc in dict(ALA='CB',ARG='CB CG CD NE CZ NH1 NH2',ASN='CB CG OD1 ND2',ASP='CB CG OD1 OD2',CYS='CB SG',GLN='CB CG CD OE1 NE2',GLU='CB CG CD OE1 OE2',GLY='',HIS='CB CG ND1 CD2 CE1 NE2',ILE='CB CG1 CG2 CD1',LEU='CB CG CD1 CD2',LYS='CB CG CD CE NZ',MET='CB CG SD CE',PHE='CB CG CD1 CD2 CE1 CE2 CZ',PRO='CB CG CD',SER='CB OG',THR='CB OG1 CG2',TRP='CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2',TYR='CB CG CD1 CD2 CE1 CE2 CZ OH',VAL='CB CG1 CG2').items()}
REGIONS={'H1':list(range(26,36)),'H2':list(range(50,67)),'H3':list(range(97,97+H3_LENGTH)),'L1':list(range(24,35)),'L2':list(range(50,57)),'L3':list(range(89,98))}
def dump(p,obj):
 p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2)+'\n');tmp.replace(p)
def fit(x,y):
 x=np.array(x);y=np.array(y);a=x.mean(0);b=y.mean(0);u,s,vt=np.linalg.svd((x-a).T@(y-b));r=u@np.diag([1,1,np.linalg.det(u@vt)])@vt;t=b-a@r;return r,t,float(np.sqrt(np.mean(np.sum((x@r+t-y)**2,axis=1))))
def audit(m,framework_only=False):
 missing=[];bonds=[];angles=[];ss=[];gaps=[];extra_cys=[]
 for c in m:
  rr=list(c)
  for r in rr:
   absent=EXPECTED[r.resname]-{a.name for a in r}
   if absent:missing.append([c.id,r.id[1],sorted(absent)])
  for r,t in zip(rr,rr[1:]):
   if r.id[1]+1!=t.id[1]:
    gaps.append([c.id,r.id[1],t.id[1]]);continue
   d=float(r['C']-t['N']);an1=float(np.degrees(calc_angle(r['CA'].get_vector(),r['C'].get_vector(),t['N'].get_vector())));an2=float(np.degrees(calc_angle(r['C'].get_vector(),t['N'].get_vector(),t['CA'].get_vector())));om=float(np.degrees(calc_dihedral(r['CA'].get_vector(),r['C'].get_vector(),t['N'].get_vector(),t['CA'].get_vector())))
   if not 1.1<=d<=1.6:bonds.append([c.id,r.id[1],d])
   if not 100<=an1<=140 or not 100<=an2<=140 or abs(om)<140 and t.resname!='PRO':angles.append([c.id,r.id[1],an1,an2,om])
  if c.id in 'HL':
   # This campaign uses the audited hu4D5 motif: validate its conserved pair,
   # independently of provisional CDR identities that AbMPNN will redesign.
   positions={'H':(22,96),'L':(23,88)}[c.id]
   sg=[c[n]['SG'].coord for n in positions if n in c and c[n].resname=='CYS' and 'SG' in c[n]]
   ss.append([c.id,float(np.linalg.norm(sg[0]-sg[1])) if len(sg)==2 else None])
   extra_cys.extend([c.id,r.id[1]] for r in c if r.resname=='CYS' and r.id[1] not in positions)
 aa=[(c.id,r.id[1],a) for c in m for r in c for a in r if a.element not in ['H','D']];clashes=[]
 for i,j in cKDTree([a.coord for c,n,a in aa]).query_pairs(1.8):
  c,n,a=aa[i];d,k,b=aa[j]
  if c==d and abs(n-k)<=1 or a.name==b.name=='SG':continue
  clashes.append([c,n,a.name,d,k,b.name,float(a-b)])
 allowed_gaps=framework_only and all(c in 'HL' for c,n,k in gaps)
 return dict(audit_version='canonical_framework_disulfides_v2',extra_CDR_cysteines=extra_cys,missing_atoms=missing,peptide_bonds=bonds,peptide_angles=angles,gaps=gaps,disulfides=ss,internal_clashes_under1p8=clashes,pass_gate=not missing and not bonds and not angles and not clashes and (not gaps or allowed_gaps) and all(x[1] is not None and 1.9<x[1]<2.2 for x in ss))
