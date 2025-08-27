//  Copyright (C) 2025 Eisuke Kawashima
//
//  Chemical JSON
//  https://github.com/OpenChemistry/chemicaljson
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
#include <string>
#include <boost/json/src.hpp>
// #include <RDGeneral/Invariant.h>

namespace RDKit {

void MolToCJSONBlock(std::ostream &os, const ROMol &mol,
                     const CJSONWriterParams &params, int confId) {
  constexpr unsigned chemicalJson = 1;  // version
  const Conformer *conf =
      mol.getNumConformers() ? &mol.getConformer(confId) : nullptr;
  const auto is3D = conf != nullptr && conf->is3D();

  const auto nAtoms = mol.getNumAtoms();
  int totalCharge = 0;

  boost::json::array coords_3d;
  boost::json::array elements_number;
  boost::json::array formalCharges;
  boost::json::array bond_indices, bond_orders;
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
      switch (params.coords) {
        case CJSONCoords::_3d:
          coords_3d.push_back(pos.x);
          coords_3d.push_back(pos.y);
          coords_3d.push_back(pos.z);
          break;
        case CJSONCoords::_3dFractional:
          throw std::invalid_argument{
              "atoms.coords.3dFractional is not supported yet"};
      }
    }
  }

  RWMol rwmol{mol};
  // boost::dynamic_bitset<> aromaticBonds;
  if (params.kekulize) {
    for (const auto bond : rwmol.bonds()) {
      if (bond->getIsAromatic()) {
        // aromaticBonds.set(bond->getIdx());
      }
      MolOps::Kekulize(rwmol);
    }
  }

  for (auto atom_itr = rwmol.beginAtoms(), atom_itr_end = rwmol.endAtoms();
       atom_itr != atom_itr_end; ++atom_itr) {
    const auto &atom = *atom_itr;
    const auto src = atom->getIdx();
    for (auto bond_itrs = rwmol.getAtomBonds(atom);
         bond_itrs.first != bond_itrs.second; ++bond_itrs.first) {
      auto *bptr = rwmol[*bond_itrs.first];
      auto *nptr = bptr->getOtherAtom(atom);
      const auto dst = nptr->getIdx();
      if (dst < src) {
        continue;
      }
      bond_indices.push_back(src);
      bond_indices.push_back(dst);

      const auto btype = bptr->getBondType();
      switch (btype) {
        case Bond::ZERO:
          bond_orders.push_back(0);
          break;
        case Bond::SINGLE:
          bond_orders.push_back(1);
          break;
        case Bond::DOUBLE:
          bond_orders.push_back(2);
          break;
        case Bond::TRIPLE:
          bond_orders.push_back(3);
          break;
        case Bond::QUADRUPLE:
          bond_orders.push_back(4);
          break;
        case Bond::QUINTUPLE:
          bond_orders.push_back(5);
          break;
        case Bond::HEXTUPLE:
          bond_orders.push_back(6);
          break;

        case Bond::DATIVE:
          [[fallthrough]];
        case Bond::DATIVEL:
          [[fallthrough]];
        case Bond::DATIVER:
          bond_orders.push_back(1);
          break;

        // XXX MUST BE INTEGER
        case Bond::AROMATIC:
          bond_orders.push_back(1.5);
          break;

        case Bond::DATIVEONE:
          bond_orders.push_back(0.5);
          break;

        case Bond::THREECENTER:
          bond_orders.push_back(0.5);
          break;
        case Bond::ONEANDAHALF:
          bond_orders.push_back(1.5);
          break;
        case Bond::TWOANDAHALF:
          bond_orders.push_back(2.5);
          break;
        case Bond::THREEANDAHALF:
          bond_orders.push_back(3.5);
          break;
        case Bond::FOURANDAHALF:
          bond_orders.push_back(4.5);
          break;
        case Bond::FIVEANDAHALF:
          bond_orders.push_back(5.5);
          break;

          // default:
          //   BOOST_LOG(rdInfoLog)
          //       << boost::format{"CMLWriter: Unsupported BondType %1%\n"} %
          //       btype;
          //   bond.put("<xmlattr>.order", "unknown");
      }
    }
  }

  std::string name;
  mol.getPropIfPresent(common_properties::_Name, name);

  boost::json::object atoms;
  atoms["coords"].emplace_object()["3d"] = coords_3d;
  atoms["elements"].emplace_object()["number"] = elements_number;
  atoms["formalCharges"] = formalCharges;

  boost::json::object bonds;
  bonds["connections"].emplace_object()["index"] = bond_indices;
  bonds["order"] = bond_orders;

  boost::json::object properties;
  properties["totalCharge"] = totalCharge;
  // properties["totalSpinMultiplicity"] = 1u;

  boost::json::object root;
  root["chemicalJson"] = chemicalJson;
  if (!name.empty()) {
    root["name"] = name;
  }
  root["atoms"] = atoms;
  root["bonds"] = bonds;
  root["properties"] = properties;
  os << root << '\n';
}

std::string MolToCJSONBlock(const ROMol &mol, const CJSONWriterParams &params,
                            int confId) {
  std::ostringstream ss;
  MolToCJSONBlock(ss, mol, params, confId);
  return ss.str();
}

void MolToCJSONFile(const ROMol &mol, const std::string &fName,
                    const CJSONWriterParams &params, int confId) {
  std::ofstream ofs{fName};
  MolToCJSONBlock(ofs, mol, params, confId);
}
}  // namespace RDKit
