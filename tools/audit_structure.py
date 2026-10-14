"""Full-heavy-atom geometry audit with campaign-specific antibody numbering."""
import argparse
from pathlib import Path
from Bio.PDB import PDBParser
from geometry import audit,dump
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 m=PDBParser(QUIET=True).get_structure('input',a.input)[0];result=audit(m);dump(a.output,result)
 print('PASS' if result['pass_gate'] else 'FAIL');raise SystemExit(0 if result['pass_gate'] else 1)
