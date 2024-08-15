# $Id$
#
#  Copyright (C) 2006  Greg Landrum
#
import os
import pickle
import sys
import unittest

from rdkit import Chem, DataStructs, RDConfig
from rdkit.Chem import MolCatalog


class TestCase(unittest.TestCase):

  def test1(self):
    cat = MolCatalog.CreateMolCatalog()
    es = []
    for smi in ('C1CCC1OC', 'C1CCC1', 'C'):
      m = Chem.MolFromSmiles(smi)
      entry = MolCatalog.MolCatalogEntry()
      entry.SetMol(m)
      self.assertTrue(entry.GetMol())
      eSmi = Chem.MolToSmiles(entry.GetMol())
      self.assertEqual(eSmi, Chem.MolToSmiles(m))
      entry.SetDescription(smi)
      self.assertEqual(entry.GetDescription(), smi)
      es.append(entry)

    v = cat.AddEntry(es[0])
    self.assertEqual(v, 0)
    self.assertEqual(cat.GetNumEntries(), 1)

    v = cat.AddEntry(es[1])
    self.assertEqual(v, 1)
    self.assertEqual(cat.GetNumEntries(), 2)

    v = cat.AddEntry(es[2])
    self.assertEqual(v, 2)
    self.assertEqual(cat.GetNumEntries(), 3)

    cat.AddEdge(0, 1)
    cat.AddEdge(0, 2)
    cat.AddEdge(1, 2)

    d = pickle.dumps(cat)
    es = None
    entry = None
    cat = None

    cat = pickle.loads(d)
    self.assertEqual(cat.GetNumEntries(), 3)
    cat = None


if __name__ == '__main__':
  unittest.main()
