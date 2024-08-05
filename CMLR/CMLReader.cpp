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

  // http://www.xml-cml.org/convention/molecular#molecule-atomArray
  // > A molecule MAY contain a single atomArray child except when it contains
  // > child molecules.
  if (auto aa = get_array("atomArray")) {
    parse_atomArray(*aa);
  }

  // http://www.xml-cml.org/convention/molecular#molecule-bondArray
  // > A molecule MAY contain a single bondyArray child provided that it does
  // > not contain child molecules.
  if (auto ba = get_array("bondArray")) {
    parse_bondArray(*ba);
  }

  return mol;
}

void CMLMolecule::parse_atomArray(const boost::property_tree::ptree &node) {
  const auto num_atoms = static_cast<unsigned>(node.count("atom"));
  // http://www.xml-cml.org/convention/molecular#atomArray-element
  // > An atomArray element MUST contain at least one child atom element.
  if (num_atoms == 0u) {
    throw cml::CMLError{"atomArray has no atom elements"};
  }
  conformer->reserve(num_atoms);

  // unsigned atom_idx = 0u;
  for (const auto &atomitr : node) {
    if (atomitr.first == "<xmlattr>") {
      for (const auto &i : atomitr.second) {
        BOOST_LOG(rdInfoLog)
            << boost::format{R"("/@%1% (= "%2%") is ignored")"} % i.first %
                   i.second.data()
            << std::endl;
      }
      continue;
    }

    if (atomitr.first != "atom") {
      BOOST_LOG(rdInfoLog) << boost::format{"/%1% is ignored"} % atomitr.first
                           << std::endl;
      continue;
    }

    parse_atom(atomitr.second);
  }
}

void CMLMolecule::parse_atom(const boost::property_tree::ptree &node) {
  const auto id = [&]() {
    try {
      return node.get<std::string>("<xmlattr>.id");
    } catch (...) {
      throw;
    }
  }();
}

void CMLMolecule::parse_bondArray(const boost::property_tree::ptree &node) {}
void CMLMolecule::parse_bond(const boost::property_tree::ptree &node) {}
}  // namespace cml
}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit
