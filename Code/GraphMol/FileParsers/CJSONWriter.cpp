//  Copyright (C) 2025 Eisuke Kawashima
//
//  XXX Chemical JSON
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "FileParsers.h"
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <fstream>
#include <sstream>
#include <boost/json/src.hpp>
#include <RDGeneral/Invariant.h>

namespace RDKit {
void MolToCJSONBlock(std::ostream &os, const ROMol &mol, int confId,
                     bool kekulize) {
  const Conformer *conf =
      mol.getNumConformers() ? &mol.getConformer(confId) : nullptr;
  const bool is3D = conf && conf->is3D();

  const auto nAtoms = mol.getNumAtoms();

  boost::json::array coords_3d;
  boost::json::array elements_number;
  boost::json::array formalCharges;
  if (is3D) {
    coords_3d.reserve(3u * nAtoms);
  }
  elements_number.reserve(nAtoms);
  formalCharges.reserve(nAtoms);

  for (auto i = 0u; i < nAtoms; i++) {
    const auto &atom = mol.getAtomWithIdx(i);
    elements_number.push_back(atom->getAtomicNum());
    if (conf) {
      const auto &pos = conf->getAtomPos(i);
      coords_3d.push_back(pos.x);
      coords_3d.push_back(pos.y);
      coords_3d.push_back(pos.z);
    }
  }

  boost::json::object coords;
  coords["3d"].emplace_array();
  boost::json::object root;
  root["atoms"] = nullptr;
}

std::string MolToCJSONBlock(const ROMol &mol, int confId, bool kekulize) {
  std::ostringstream ss;
  MolToCJSONBlock(ss, mol, confId, kekulize);
  return ss.str();
}

void MolToCJSONFile(const ROMol &mol, const std::string &fName, int confId,
                    bool kekulize) {
  std::ofstream ofs{fName};
  MolToCJSONBlock(ofs, mol, confId, kekulize);
}
}  // namespace RDKit
