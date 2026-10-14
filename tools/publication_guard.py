"""Inspect staged publication content; never output suspected secrets or sequences."""
import argparse,re,subprocess,sys
from pathlib import PurePosixPath
ROOT_FILES={'.gitignore','README.md','DATA_POLICY.md','requirements-analysis.txt'}
SUFFIXES={'.md','.py','.sh'}
ALPHABET='ACDEFGHIKLMNPQRSTVWY'
def violations(name,data):
 p=PurePosixPath(name);issues=[]
 if not (name in ROOT_FILES or p.parts[0] in {'docs','tools','tests','hooks'} and (p.suffix in SUFFIXES or name in {'hooks/pre-commit','hooks/pre-push'})):issues.append('file not on publication allowlist')
 if len(data)>100000:issues.append('file exceeds reviewed source-size limit')
 try:text=data.decode('utf-8')
 except UnicodeDecodeError:return issues+['binary file']
 if re.search(r'(?:tnf|egfr)_[0-9a-f]{16}',text):issues.append('private candidate identifier')
 if any(m.group()!=ALPHABET for m in re.finditer(r'[ACDEFGHIKLMNPQRSTVWY]{12,}',text)):issues.append('possible protein sequence literal')
 if re.search(r'^(?:ATOM  |HETATM)\s*\d',text,re.M):issues.append('structural coordinate record')
 if re.search(r'(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,}|-----BEGIN [A-Z ]*PRIVATE KEY-----)',text):issues.append('credential-like material')
 if re.search(r'/(?:data[0-9]+|home)/[a-z][a-z0-9_-]+/',text):issues.append('private server path')
 return issues

def run(ref=None):
 cmd=['git','ls-tree','-r','--name-only','-z',ref] if ref else ['git','ls-files','-z']
 names=subprocess.check_output(cmd).decode().split('\0');bad=[]
 for name in filter(None,names):
  mode=subprocess.check_output(['git','ls-tree',ref,'--',name]).decode().split()[0] if ref else subprocess.check_output(['git','ls-files','--stage','--',name]).decode().split()[0]
  issues=['symlink/submodule prohibited'] if mode in {'120000','160000'} else []
  data=subprocess.check_output(['git','show',f'{ref}:{name}' if ref else ':'+name]);issues+=violations(name,data)
  if issues:bad.append((name,issues))
 for name,issues in bad:print(name+': '+', '.join(issues),file=sys.stderr)
 print(f'Publication guard: {len(list(filter(None,names)))} files; {len(bad)} flagged')
 return bool(bad)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--ref');a=p.parse_args();raise SystemExit(run(a.ref))
