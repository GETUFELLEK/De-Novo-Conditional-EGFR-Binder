"""Descriptive sequence motifs; no single motif establishes a liability or rejection."""
import argparse,csv,re
from pathlib import Path
MOTIFS={'N_linked_sequon':r'N[^P][ST]','NG_or_NS':r'N[GS]','DG_or_DP':r'D[GP]','hydrophobic_run_5':r'[AVILMFWY]{5,}','glycine_run_5':r'G{5,}'}
def flags(sequence):
 if not sequence or re.search(r'[^ACDEFGHIKLMNPQRSTVWY]',sequence):raise ValueError('Expected a nonempty canonical protein sequence')
 return [(label,m.start()+1,m.group(1)) for label,pattern in MOTIFS.items() for m in re.finditer('(?=('+pattern+'))',sequence)]
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 with a.input.open() as f:rows=list(csv.DictReader(f))
 a.output.parent.mkdir(parents=True,exist_ok=True)
 with a.output.open('w') as f:
  w=csv.writer(f);w.writerow(['name','chain','flag','chain_local_start','motif'])
  for x in rows:
   for chain,column in [('H','VH'),('L','VL')]:
    for label,pos,motif in flags(x[column]):w.writerow([x['name'],chain,label,pos,motif])
