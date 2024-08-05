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

  unsigned atom_idx = 0u;
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

    parse_atom(atomitr.second, atom_idx++);
  }
}

void CMLMolecule::parse_atom(const boost::property_tree::ptree &node,
                             unsigned idx) {
  const auto id = [&]() {
    std::string tmp;
    // http://www.xml-cml.org/convention/molecular#atom-id
    // > An atom MUST have an id attribute
    try {
      tmp = node.get<std::string>("<xmlattr>.id");
    } catch (const boost::property_tree::ptree_bad_path &) {
      throw cml::CMLError{"atom/@id is missing"};
    }
    // http://www.xml-cml.org/convention/molecular#atom-id
    // > The value of the id MUST be unique amongst the atoms within the eldest
    // > containing molecule.
    if (id_atom.count(tmp)) {
      auto msg = boost::format{R"("atom/@id (= "%1%") is not unique")"} % tmp;
      throw cml::CMLError{msg.str()};
    }
    return tmp;
  }();
  id_atom[id] = nullptr;

  auto atom = [&]() -> std::unique_ptr<RDKit::Atom> {
    std::string tmp;
    // http://www.xml-cml.org/convention/molecular#atom-elementType
    // > An atom MUST have an elementType attribute
    try {
      tmp = node.get<std::string>("<xmlattr>.elementType");
    } catch (const boost::property_tree::ptree_bad_path &) {
      throw cml::CMLError{"atom/@elementType is missing"};
    }

    // http://www.xml-cml.org/schema/schema3/schema.xsd
    // > Du. ("dummy") This does not correspond to a "real" atom and can support
    // > a point in space or within a chemical graph. R. ("R-group") This
    // > indicates that an atom or group of atoms could be > attached at this
    // > point.
    return (tmp == "Du" || tmp == "R") ? std::make_unique<RDKit::Atom>()
                                       : std::make_unique<RDKit::Atom>(tmp);
  }();

  // http://www.xml-cml.org/convention/molecular#atom-isotopeNumber
  // > An atom MAY have an isotopeNumber attribute.
  if (auto isotopeNumber =
          node.get_optional<unsigned>("<xmlattr>.isotopeNumber")) {
    atom->setIsotope(*isotopeNumber);
  }

  // http://www.xml-cml.org/convention/molecular#atom-formalCharge
  // > An atom MAY have a formalCharge attribute.
  if (auto formalCharge = node.get_optional<int>("<xmlattr>.formalCharge")) {
    atom->setFormalCharge(*formalCharge);
    sum_formalCharge += *formalCharge;
  }

  // http://www.xml-cml.org/convention/molecular#atom-spinMultiplicity
  // > An atom MAY have a spinMultiplicity attribute.
  if (auto spinMultiplicity =
          node.get_optional<unsigned>("<xmlattr>.spinMultiplicity")) {
    if (*spinMultiplicity == 0u) {
      auto msg = boost::format{"atom/@spinMultiplicity is zero"};
      throw cml::CMLError{msg.str()};
    }
    if (*spinMultiplicity <= 2u) {
      atom->setNumRadicalElectrons(*spinMultiplicity - 1u);
    } else {
      // TODO
      auto msg = boost::format{"atom/@spinMultiplicity (= %2%) is ignored"} %
                 *spinMultiplicity;
      BOOST_LOG(rdWarningLog) << msg;
    }
  }

  const auto hydrogenCount =
      node.get_optional<unsigned>("<xmlattr>.hydrogenCount");
  id_hydrogenCount[id] = hydrogenCount;

  // http://www.xml-cml.org/convention/molecular#atom-x3
  // > An atom MAY have an x3 attribute, the value of which is the x coordinate
  // > of a 3 dimensional object. The units are Angstrom and the axis system is
  // > always right handed.
  const auto x3 = node.get_optional<double>("<xmlattr>.x3");
  const auto y3 = node.get_optional<double>("<xmlattr>.y3");
  const auto z3 = node.get_optional<double>("<xmlattr>.z3");

  if (x3 || y3 || z3) {
    // http://www.xml-cml.org/convention/molecular#atom-x3
    // > If a x3 attribute is present there MUST also be a y3 and z3 present.
    if (!(x3 && y3 && z3)) {
      auto msg =
          boost::format{"/atom/ does not have all of x3, y3 and z3 attributes"};
      throw cml::CMLError{msg.str()};
    }
    RDGeom::Point3D r{*x3, *y3, *z3};
    // BOOST_LOG(rdDebugLog) << xpath_to_atom << ' ' << r << std::endl;
    conformer->setAtomPos(idx, r);
  } else {
    // Ignore x2 and y2
    // http://www.xml-cml.org/convention/molecular#atom-x2
    // > An atom MAY have an x2 attribute, the value of which is used for
    // > displaying the object in 2 dimensions. This is unrelated to the 3-D
    // > coordinates for the object.
    BOOST_LOG(rdInfoLog)
        << boost::format{"/atom does not have geometrical info "
                 "(x2 and y2 are ignored if exist)"}
        << std::endl;
  }

  id_atom[id] = std::move(atom);
}

void CMLMolecule::parse_bondArray(const boost::property_tree::ptree &node) {}
void CMLMolecule::parse_bond(const boost::property_tree::ptree &node) {}
}  // namespace cml
}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit
