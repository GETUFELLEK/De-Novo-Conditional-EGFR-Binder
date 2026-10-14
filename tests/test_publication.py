"""Publication-boundary tests use constructed tokens, never real design data."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from publication_guard import violations
from sequence_flags import flags
class BoundaryTests(unittest.TestCase):
 def test_safe_document(self):self.assertEqual(violations('docs/methods.md',b'# Methods\nRead-only description.'),[])
 def test_rejects_coordinate_files(self):self.assertTrue(violations('tools/input.pdb',b'content'))
 def test_rejects_embedded_coordinates(self):self.assertTrue(violations('docs/test.md',b'ATOM      1  CA'))
 def test_rejects_sequence_literal(self):self.assertTrue(violations('tools/example.py',('A'*30).encode()))
 def test_rejects_candidate_identifier(self):self.assertTrue(violations('docs/test.md',('tnf_'+'0'*16).encode()))
 def test_rejects_token(self):self.assertTrue(violations('docs/test.md',('ghp_'+'0'*30).encode()))
 def test_sequon_excludes_proline(self):
  self.assertIn(('N_linked_sequon',1,'NAT'),flags('NAT'))
  self.assertNotIn('N_linked_sequon',[x[0] for x in flags('NPT')])
 def test_invalid_sequence(self):
  with self.assertRaises(ValueError):flags('X')
if __name__=='__main__':unittest.main()
