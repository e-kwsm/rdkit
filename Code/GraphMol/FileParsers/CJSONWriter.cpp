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
#include <boost/format.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <RDGeneral/Invariant.h>

namespace RDKit {
void MolToCJSONBlock(std::ostream &os, const ROMol &mol, int confId,
                   bool kekulize) {
  auto pt = molToPTree(mol, confId, kekulize);
  if (pt.empty()) {
    return;
  }
  boost::property_tree::write_xml(
      os, pt,
      boost::property_tree::xml_writer_make_settings<std::string>(' ', 2));
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
}//
