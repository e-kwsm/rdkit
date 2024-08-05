//  Copyright (C) 2020 Eisuke Kawashima
//
//  Chemical Markup Language, CML, schema 3
//  http://www.xml-cml.org/
//  See
//  http://www.xml-cml.org/convention/molecular
//  http://www.xml-cml.org/schema/schema3/schema.xsd
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "CMLReader.h"
#include <GraphMol/RWMol.h>

namespace RDKit {
namespace v2 {
namespace FileParsers {
namespace cml {
CMLError::~CMLError() = default;

CMLMolecule::CMLMolecule(const boost::property_tree::ptree &molecule_node)
    : molecule_node{molecule_node},
      molecule{std::make_unique<RDKit::RWMol>()},
      conformer{std::make_unique<RDKit::Conformer>()} {}

CMLMolecule::~CMLMolecule() = default;

std::unique_ptr<RWMol> CMLMolecule::parse() {
  auto mol = std::make_unique<RWMol>();
  return mol;
}

void CMLMolecule::parse_atomArray(const boost::property_tree::ptree &node) {}
void CMLMolecule::parse_atom(const boost::property_tree::ptree &node) {}
void CMLMolecule::parse_bondArray(const boost::property_tree::ptree &node) {}
void CMLMolecule::parse_bond(const boost::property_tree::ptree &node) {}
}  // namespace cml
}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit
