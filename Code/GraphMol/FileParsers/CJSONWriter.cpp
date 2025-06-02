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

#include "CJson.h"
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
  const auto is3D = conf != nullptr && conf->is3D();

  const auto nAtoms = mol.getNumAtoms();
  int totalCharge = 0;

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
    const auto charge = atom->getFormalCharge();
    formalCharges.push_back(charge);
    totalCharge += charge;
    if (conf) {
      const auto &pos = conf->getAtomPos(i);
      coords_3d.push_back(pos.x);
      coords_3d.push_back(pos.y);
      coords_3d.push_back(pos.z);
    }
  }

  boost::json::object atoms;
  atoms["coords"].emplace_object()["3d"] = coords_3d;
  atoms["elements"].emplace_object()["number"] = elements_number;
  atoms["formalCharges"] = formalCharges;

  boost::json::object bonds;

  boost::json::object properties;
  properties["totalCharge"] = totalCharge;
  properties["totalSpinMultiplicity"] = 1u;

  boost::json::object root;
  root["atoms"] = atoms;
  root["bonds"] = bonds;
  root["chemicalJson"] = 1;
  root["properties"] = properties;
  os << root;
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
