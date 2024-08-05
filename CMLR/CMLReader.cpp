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

#include <boost/format.hpp>

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

auto CMLMolecule::get_array(const std::string &name) const {
  std::unique_ptr<boost::property_tree::ptree> p;
  switch (molecule_node.count(name)) {
    case 0u:
      BOOST_LOG(rdInfoLog) << boost::format{"%s is missing"} % name
                           << std::endl;
      break;
    case 1u:
      p = std::make_unique<boost::property_tree::ptree>(
          molecule_node.get_child(name));
      break;
    default:
      auto msg = boost::format{"multiple %s elements"} % name;
      throw cml::CMLError{msg.str()};
  }
  return p;
}

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
