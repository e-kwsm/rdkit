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

using namespace std::literals::string_literals;

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
      throw CMLError{msg.str()};
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
    throw CMLError{"atomArray has no atom elements"};
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
      throw CMLError{"atom/@id is missing"};
    }
    // http://www.xml-cml.org/convention/molecular#atom-id
    // > The value of the id MUST be unique amongst the atoms within the eldest
    // > containing molecule.
    if (id_atom.count(tmp)) {
      auto msg = boost::format{R"("atom/@id (= "%1%") is not unique")"} % tmp;
      throw CMLError{msg.str()};
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
      throw CMLError{"atom/@elementType is missing"};
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
      throw CMLError{msg.str()};
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
      throw CMLError{msg.str()};
    }
    RDGeom::Point3D r{*x3, *y3, *z3};
    conformer->setAtomPos(idx, r);
  } else {
    // Ignore x2 and y2
    // http://www.xml-cml.org/convention/molecular#atom-x2
    // > An atom MAY have an x2 attribute, the value of which is used for
    // > displaying the object in 2 dimensions. This is unrelated to the 3-D
    // > coordinates for the object.
    BOOST_LOG(rdInfoLog)
        << boost::format{"/atom does not have geometrical info "
                 "(x2 and y2 are ignored even if exist)"}
        << std::endl;
  }

  id_atom[id] = std::move(atom);
}

void CMLMolecule::parse_bondArray(const boost::property_tree::ptree &node) {
  // http://www.xml-cml.org/convention/molecular#bondArray-element
  // > A bondArray element MUST contain at least one child bond element.
  if (node.count("bond") == 0u) {
    auto msg = boost::format{"/bondArray has no bond elements"};
    throw CMLError{msg.str()};
  }

  for (const auto &bonditr : node) {
    if (bonditr.first == "<xmlattr>") {
      continue;
    }

    if (bonditr.first != "bond") {
      BOOST_LOG(rdInfoLog) << boost::format{"bondArray/%1% is ignored"} %
                                  bonditr.first
                           << std::endl;
      continue;
    }

    parse_bond(bonditr.second);
  }
}

void CMLMolecule::parse_bond(const boost::property_tree::ptree &node) {
  auto xpath_to_bond = "/bond"s;
  // http://www.xml-cml.org/convention/molecular#bond-id
  // > It is RECOMMENDED that a bond has an id attribute so that it can be
  // > referenced.
  const auto id = node.get_optional<std::string>("<xmlattr>.id");
  if (!id) {
    BOOST_LOG(rdInfoLog) << boost::format{"%1%/@id is missing"} % xpath_to_bond
                         << std::endl;
  }
  // } else if (!is_valid_id(*id)) {
  //   auto msg =
  //       boost::format{"%1%/@id (= \"%2%\") is invalid"} % xpath_to_bond %
  //       id;
  //   throw RDKit::FileParseException{msg.str()};
  else if (!bond_ids.insert(*id).second) {
    // http://www.xml-cml.org/convention/molecular#bond-id
    // > The id of a bond MUST be unique amongst the bonds of the eldest
    // > containing molecule.
    auto msg = boost::format{"%1%/@id (= \"%2%\") is not unique"} %
               xpath_to_bond % *id;
    throw RDKit::FileParseException{msg.str()};
  } else {
    xpath_to_bond += "[@id=\"" + *id + "\"]";
  }

  // http://www.xml-cml.org/convention/molecular#bond-order
  // > A bond MUST have an order attribute
  const auto order = [&]() {
    std::string tmp;
    try {
      tmp = node.get<std::string>("<xmlattr>.order");
    } catch (...) {
      auto msg = boost::format{"%1%/@order is missing"} % xpath_to_bond;
      throw CMLError{msg.str()};
    }
    return tmp;
  }();

  RDKit::Bond::BondType bt = RDKit::Bond::UNSPECIFIED;
  if (order == "1" || order == "S") {
    bt = RDKit::Bond::SINGLE;
  } else if (order == "2" || order == "D") {
    bt = RDKit::Bond::DOUBLE;
  } else if (order == "3" || order == "T") {
    bt = RDKit::Bond::TRIPLE;
  } else if (order == "A") {
    bt = RDKit::Bond::AROMATIC;
  } /* RDKit extension */ else if (order == "4") {
    bt = RDKit::Bond::QUADRUPLE;
  } else if (order == "5") {
    bt = RDKit::Bond::QUINTUPLE;
  } else if (order == "6") {
    bt = RDKit::Bond::HEXTUPLE;
  } else if (order == "1.5") {
    bt = RDKit::Bond::ONEANDAHALF;
  } else if (order == "2.5") {
    bt = RDKit::Bond::TWOANDAHALF;
  } else if (order == "3.5") {
    bt = RDKit::Bond::THREEANDAHALF;
  } else if (order == "4.5") {
    bt = RDKit::Bond::FOURANDAHALF;
  } else if (order == "5.5") {
    bt = RDKit::Bond::FIVEANDAHALF;
  } else {
    BOOST_LOG(rdWarningLog)
        << boost::format{R"("%1%/@order (= "%2%") is unrecognizable")"} %
               xpath_to_bond % order
        << std::endl;
  }

  // http://www.xml-cml.org/convention/molecular#bond-atomRefs2
  // > A bond MUST have a atomRefs2 attribute, the value of which MUST be the
  // > space separated ids of two different atoms which MUST be in the same
  // > molecule.
  const auto atomRefs2 = [&]() {
    std::string tmp;
    try {
      tmp = node.get<std::string>("<xmlattr>.atomRefs2");
    } catch (...) {
      auto msg = boost::format{"%1%/@atomRefs2 is missing"} % xpath_to_bond;
      throw CMLError{msg.str()};
    }
    return tmp;
  }();

  std::istringstream iss{atomRefs2};
  std::string id_bgn, id_end;
  if (!(iss >> id_bgn >> id_end)) {
    auto msg =
        boost::format{R"(%1%/@atomRefs2 (= "%2%") does not have two ids ")"
                      "separated by space"} %
        xpath_to_bond % atomRefs2;
    throw CMLError{msg.str()};
  } else if (std::string extra; iss >> extra) {
    auto msg =
        boost::format{R"(%1%/@atomRefs2 (= "%2%") has three or more ids)"} %
        xpath_to_bond % atomRefs2;
    throw CMLError{msg.str()};
  }

  if (id_bgn == id_end) {
    auto msg = boost::format{R"(%1%/@atomRefs2 (= "%2%") is self-bond)"} %
               xpath_to_bond % atomRefs2;
    throw CMLError{msg.str()};
  }

  auto distance = [&](const auto &atom_id) {
    auto i = id_atom.find(atom_id);
    if (i == id_atom.end()) {
      auto msg =
          boost::format{
              R"(%1%/@atomRefs2 (= \"%2%\") refers to non-existing )"
              R"(../../../atomArray/atom[@id="%3%"])"} %
          xpath_to_bond % atomRefs2 % atom_id;
      throw CMLError{msg.str()};
    }
    return static_cast<unsigned>(std::distance(id_atom.begin(), i));
  };

  const auto index_bgn = distance(id_bgn);
  const auto index_end = distance(id_end);

  auto bond = std::make_unique<RDKit::Bond>(bt);
  bond->setBeginAtomIdx(index_bgn);
  bond->setEndAtomIdx(index_end);

  // TODO
  const auto bondStereo = node.get_optional<std::string>("bondStereo");
  if (bondStereo) {
    if (*bondStereo == "C") {
      // TODO
      // b->setStereo(Bond::STEREOCIS);
    } else if (*bondStereo == "T") {
      // TODO
      // b->setStereo(Bond::STEREOTRANS);
    } else if (*bondStereo == "W") {
    } else if (*bondStereo == "H") {
    } else if (*bondStereo == "other") {
      BOOST_LOG(rdWarningLog) << "not implemented" << std::endl;
    } else if (*bondStereo != "undefined") {
      BOOST_LOG(rdWarningLog)
          << boost::format{"%1%/bondStereo (= \"%2%\") is unrecognizable"} %
                 xpath_to_bond % *bondStereo
          << std::endl;
    }
  }

  bond->setOwningMol(*molecule);
  molecule->addBond(bond.release(), true);
}
}  // namespace cml
}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit
