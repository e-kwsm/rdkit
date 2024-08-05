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

// TODO namespace support

#include <optional>
#ifndef __clang__version__                     // FIXME
#pragma GCC diagnostic push                    // FIXME
#pragma GCC diagnostic ignored "-Wall"         // FIXME
#pragma GCC diagnostic ignored "-Wextra"       // FIXME
#pragma GCC diagnostic ignored "-Wpedantic"    // FIXME
#else                                          // FIXME
#pragma clang diagnostic push                  // FIXME
#pragma clang diagnostic ignored "-Wall"       // FIXME
#pragma clang diagnostic ignored "-Wextra"     // FIXME
#pragma clang diagnostic ignored "-Wpedantic"  // FIXME
#endif                                         // FIXME

#include <fstream>
#include <iostream>

#include <boost/exception/diagnostic_information.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/xml_parser.hpp>

// #include "FileParsers.h"
#include "CMLReader.h"
#include <GraphMol/RWMol.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/BadFileException.h>

#ifndef __clang__version__    // FIXME
#pragma GCC diagnostic pop    // FIXME
#else                         // FIXME
#pragma clang diagnostic pop  // FIXME
#endif                        // FIXME

// using namespace std::literals::string_literals;  // XXX

namespace RDKit {
namespace v2 {
namespace FileParsers {
namespace cml {
CMLError::~CMLError() = default;
XMLMalformedError::~XMLMalformedError() = default;
MandatoryElementNotFoundError::~MandatoryElementNotFoundError() = default;
CMLAttributeError::~CMLAttributeError() = default;

CMLMolecule::CMLMolecule(const boost::property_tree::ptree &molecule_node)
    : molecule_node{molecule_node},
      molecule{std::make_unique<RDKit::RWMol>()},
      conformer{std::make_unique<RDKit::Conformer>()} {}

CMLMolecule::~CMLMolecule() = default;

std::unique_ptr<RWMol> CMLMolecule::parse() {
  auto check_molecule = []() {};          // XXX
  auto xpath = [](auto x) { return x; };  // XXX
  const auto molecule_xpath = __func__;

  check_molecule();

  // http://www.xml-cml.org/convention/molecular#molecule-name
  // > A molecule MAY contain any number of name children.
  // XXX take the first occurrence
  if (auto name = molecule_node.get_optional<std::string>("name")) {
    molecule->setProp(common_properties::_Name, *name);
  }

  auto get_array = [&](auto node) {
    std::unique_ptr<boost::property_tree::ptree> p;
    switch (molecule_node.count(node)) {
      case 0u:
        BOOST_LOG(rdInfoLog)
            << boost::format{"%1% is missing"} % xpath(node) << std::endl;
        break;
      case 1u:
        p = std::make_unique<boost::property_tree::ptree>(
            molecule_node.get_child(node));
        break;
      default:
        auto msg = boost::format{"%1% has multiple %2% elements"} %
                   molecule_xpath % node;
        throw RDKit::FileParseException{msg.str()};
    }
    return p;
  };

  // http://www.xml-cml.org/convention/molecular#molecule-atomArray
  // > A molecule MAY contain a single atomArray child except when it contains
  // > child molecules.
  if (auto aa = get_array("atomArray")) {
    // parse_atomArray(xpath("atomArray").str(), *aa);
  }

  // http://www.xml-cml.org/convention/molecular#molecule-bondArray
  // > A molecule MAY contain a single bondyArray child provided that it does
  // > not contain child molecules.
  if (auto ba = get_array("bondArray")) {
    // parse_bondArray(xpath("bondArray").str(), *ba);
  }

  if (conformer) {
    molecule->addConformer(conformer.release());
  }

  //  check_hydrogenCount();
  //  if (params.sanitize) {
  //    if (params.removeHs) {
  //      MolOps::RemoveHsParameters ps;
  //      MolOps::removeHs(*molecule, ps, true);
  //    } else {
  //      MolOps::sanitizeMol(*molecule);
  //    }
  //  } else {
  //    // molecule->updatePropertyCache(false);
  //  }
  //  // MolOps::assignChiralTypesFrom3D(*molecule);

  return std::move(molecule);
}

void CMLMolecule::parse_atomArray(
    const boost::property_tree::ptree &atomArray) {
  const std::string xpath_to_atomArray = __func__;
  const auto num_atoms = static_cast<unsigned>(atomArray.count("atom"));
  // http://www.xml-cml.org/convention/molecular#atomArray-element
  // > An atomArray element MUST contain at least one child atom element.
  if (num_atoms == 0u) {
    auto msg = boost::format{"%1% has no atom elements"} % xpath_to_atomArray;
    throw cml::MandatoryElementNotFoundError{msg.str()};
  }
  conformer->reserve(num_atoms);

  unsigned atom_idx = 0u;
  for (const auto &atomitr : atomArray) {
    if (atomitr.first == "<xmlattr>") {
      for (const auto &i : atomitr.second) {
        BOOST_LOG(rdInfoLog)
            << boost::format{"%1%/@%2% (= \"%3%\") is ignored"} %
                   xpath_to_atomArray % i.first % i.second.data()
            << std::endl;
      }
      continue;
    }

    if (atomitr.first != "atom") {
      BOOST_LOG(rdInfoLog) << boost::format{"%1%/%2% is ignored"} %
                                  xpath_to_atomArray % atomitr.first
                           << std::endl;
      continue;
    }

    parse_atom(
        // xpath_to_atomArray + "/" + atomitr.first,
        atomitr.second
        // ,atom_idx++
    );
  }

  // if (formalCharge && *formalCharge != sum_formalCharge) {
  //   auto msg = boost::format{"%1% (= %2%) is not equal to sum of %3% (=
  //   %4%)"} %
  //              xpath("@formalCharge") % *formalCharge %
  //              "../atomArray/atom/@formalCharge" % sum_formalCharge;
  //   throw RDKit::FileParseException{msg.str()};
  // }

  // TODO check spinMultiplicity

  for (auto &&p : id_atom) {
    molecule->addAtom(p.second.release(), true, true);
  }
}

void CMLMolecule::parse_atom(const boost::property_tree::ptree &atom_node) {
  std::string xpath_to_atom = __func__;
  const auto id = atom_node.get_optional<std::string>("<xmlattr>.id");
  if (!id) {
    // http://www.xml-cml.org/convention/molecular#atom-id
    // > An atom MUST have an id attribute
    auto msg = boost::format{"%1%/@id is missing"} % xpath_to_atom;
    throw RDKit::FileParseException{msg.str()};
  }
  // else if (!is_valid_id(*id)) {
  //   auto msg =
  //       boost::format{"%1%/@id (= \"%2%\") is invalid"} % xpath_to_atom %
  //       *id;
  //   throw RDKit::FileParseException{msg.str()};
  // }
  else if (id_atom.count(*id)) {
    // http://www.xml-cml.org/convention/molecular#atom-id
    // > The value of the id MUST be unique amongst the atoms within the eldest
    // > containing molecule.
    auto msg = boost::format{"%1%/@id (= \"%2%\") is not unique"} %
               xpath_to_atom % *id;
    throw RDKit::FileParseException{msg.str()};
  }
  xpath_to_atom += "[@id=\"" + *id + "\"]";

  const auto elementType =
      atom_node.get_optional<std::string>("<xmlattr>.elementType");
  // http://www.xml-cml.org/convention/molecular#atom-elementType
  // > An atom MUST have an elementType attribute
  if (!elementType) {
    auto msg = boost::format{"%1%/@elementType is missing"} % xpath_to_atom;
    throw RDKit::FileParseException{msg.str()};
  }

  // http://www.xml-cml.org/schema/schema3/schema.xsd
  // > Du. ("dummy") This does not correspond to a "real" atom and can support a
  // > point in space or within a chemical graph.
  // > R. ("R-group") This indicates that an atom or group of atoms could be
  // > attached at this point.
  auto atom = (*elementType == "Du" || *elementType == "R")
                  ? std::make_unique<RDKit::Atom>()
                  : std::make_unique<RDKit::Atom>(*elementType);

#if 0
  // http://www.xml-cml.org/convention/molecular#atom-isotopeNumber
  // > An atom MAY have an isotopeNumber attribute.
  if (auto isotopeNumber = atom_node.get_optional<unsigned>(
          "<xmlattr>.isotopeNumber",
          LexicalTranslator<unsigned>{xpath_to_atom + "/@isotopeNumber"})) {
    atom->setIsotope(*isotopeNumber);
  }

  // http://www.xml-cml.org/convention/molecular#atom-formalCharge
  // > An atom MAY have a formalCharge attribute.
  if (auto formalCharge = atom_node.get_optional<int>(
          "<xmlattr>.formalCharge",
          LexicalTranslator<int>{xpath_to_atom + "/@formalCharge"})) {
    atom->setFormalCharge(*formalCharge);
    sum_formalCharge += *formalCharge;
  }

  // http://www.xml-cml.org/convention/molecular#atom-spinMultiplicity
  // > An atom MAY have a spinMultiplicity attribute.
  if (auto spinMultiplicity = atom_node.get_optional<unsigned>(
          "<xmlattr>.spinMultiplicity",
          LexicalTranslator<unsigned>{xpath_to_atom + "/@spinMultiplicity"})) {
    if (*spinMultiplicity == 0u) {
      auto msg = boost::format{"%1%/@spinMultiplicity is zero"} % xpath_to_atom;
      throw RDKit::FileParseException{msg.str()};
    }
    if (*spinMultiplicity <= 2u) {
      atom->setNumRadicalElectrons(*spinMultiplicity - 1u);
    } else {
      // TODO
      auto msg = boost::format{"%1%/@spinMultiplicity (= %2%) is ignored"} %
                 xpath_to_atom % *spinMultiplicity;
      BOOST_LOG(rdWarningLog) << msg;
    }
  }

  const auto hydrogenCount = atom_node.get_optional<unsigned>(
      "<xmlattr>.hydrogenCount",
      LexicalTranslator<unsigned>{xpath_to_atom + "/@hydrogenCount"});
  id_hydrogenCount[xpath_to_atom] = hydrogenCount;

  // http://www.xml-cml.org/convention/molecular#atom-x3
  // > An atom MAY have an x3 attribute, the value of which is the x coordinate
  // > of a 3 dimensional object. The units are Angstrom and the axis system is
  // > always right handed.
  const auto x3 = atom_node.get_optional<double>(
      "<xmlattr>.x3", LexicalTranslator<double>{xpath_to_atom + "/@x3"});
  const auto y3 = atom_node.get_optional<double>(
      "<xmlattr>.y3", LexicalTranslator<double>{xpath_to_atom + "/@y3"});
  const auto z3 = atom_node.get_optional<double>(
      "<xmlattr>.z3", LexicalTranslator<double>{xpath_to_atom + "/@z3"});

  if (x3 || y3 || z3) {
    // http://www.xml-cml.org/convention/molecular#atom-x3
    // > If a x3 attribute is present there MUST also be a y3 and z3 present.
    if (!(x3 && y3 && z3)) {
      auto msg =
          boost::format{"%1% does not have all of x3, y3 and z3 attributes"} %
          xpath_to_atom;
      throw RDKit::FileParseException{msg.str()};
    }
    RDGeom::Point3D r{*x3, *y3, *z3};
    BOOST_LOG(rdDebugLog) << xpath_to_atom << ' ' << r << std::endl;
    conformer->setAtomPos(idx, r);
  } else {
    // Ignore x2 and y2
    // http://www.xml-cml.org/convention/molecular#atom-x2
    // > An atom MAY have an x2 attribute, the value of which is used for
    // > displaying the object in 2 dimensions. This is unrelated to the 3-D
    // > coordinates for the object.
    BOOST_LOG(rdInfoLog)
        << boost::format{"%1% does not have geometrical info (x2 and y2 are ignored if exist)"} %
               xpath_to_atom
        << std::endl;
  }
#endif

  id_atom[*id] = std::move(atom);
}

void CMLMolecule::parse_bondArray(
    const boost::property_tree::ptree &bondArray) {}
}  // namespace cml

CMLSupplier::CMLSupplier(std::unique_ptr<std::istream> &&p_istream,
                         const CMLFileParserParams &params) noexcept(false)
    : p_istream{std::move(p_istream)}, params{params} {
  PRECONDITION(this->p_istream, "bad stream");

  try {
    boost::property_tree::read_xml(*this->p_istream, pt);
  } catch (const boost::property_tree::xml_parser_error &e) {
    throw cml::XMLMalformedError{boost::diagnostic_information(e)};
  }

  if (pt.size() > 1u) {
    throw cml::XMLMalformedError{"XML MUST NOT have multiple roots"};
  }
  const auto root = pt.front();
  if (root.first == "molecule") {
    molecule_node = root.second;
  } else {
    // http://www.xml-cml.org/convention/molecular#applying
    // > If the molecular convention is specified on a cml element then that
    // > element MUST have at least one child molecule element that either has
    // > no convention specified or specifies the molecular convention.
    if (root.first != "cml") {
      BOOST_LOG(rdWarningLog)
          << boost::format{"XML root element is %s"} % root.first << std::endl;
    }
    auto m = root.second.find("molecule");
    // XXX do not dig into 3rd or deeper
    if (m == root.second.not_found()) {
      BOOST_LOG(rdWarningLog)
          << boost::format{"/%s/molecule is not found"} % root.first
          << std::endl;
    } else {
      molecule_node = m->second;
    }
  }
}

CMLSupplier::CMLSupplier(const std::string &fileName,
                         const CMLFileParserParams &params) noexcept(false)
    : params{params} {
  p_istream = std::make_unique<std::ifstream>(fileName);
  if (!p_istream || !p_istream->good()) {
    auto msg = boost::format{"Bad input file %1%"} % fileName;
    throw BadFileException{msg.str()};
  }
}

CMLSupplier::~CMLSupplier() = default;

void CMLSupplier::close() { p_istream.reset(); }

std::unique_ptr<RWMol> CMLSupplier::next() {
  if (!molecule_node) {
    throw cml::CMLError{"EOF"};
  }
  // auto tmp = parse_molecule_node(molecule_node.value());
  auto tmp = cml::CMLMolecule{molecule_node.value()}.parse();
  molecule_node = std::nullopt;  // XXX
  return tmp;
}

std::unique_ptr<boost::property_tree::ptree> CMLSupplier::get_array(
    const boost::property_tree::ptree &node, const std::string &name) const
    noexcept(false) {
  auto xpath = [](auto x) { return x; };
  switch (node.count(name)) {
    case 0u:
      BOOST_LOG(rdInfoLog) << boost::format{"%1% is missing"} % xpath(name)
                           << std::endl;
      return nullptr;
    case 1u:
      return std::make_unique<boost::property_tree::ptree>(
          node.get_child(name));
    default:
      // auto msg = boost::format{"%1% has multiple %2% elements"} %
      // molecule_xpath % node; throw RDKit::FileParseException{msg.str()};
      throw;
  }
}

std::unique_ptr<RWMol> CMLSupplier::parse_molecule_node(
    const boost::property_tree::ptree &node) {
  std::unique_ptr<RWMol> mol;

  constexpr auto xpath_to_atomArray = __func__;
  constexpr auto xpath_to_atom = __func__;

  // http://www.xml-cml.org/convention/molecular#molecule-atomArray
  // > A molecule MAY contain a single atomArray child except when it contains
  // > child molecules.
  if (auto atomArray = get_array(node, "atomArray")) {
    // parse_atomArray(xpath("atomArray").str(), *aa);
    const auto num_atoms = static_cast<unsigned>(atomArray->count("atom"));
    // http://www.xml-cml.org/convention/molecular#atomArray-element
    // > An atomArray element MUST contain at least one child atom element.
    if (num_atoms == 0u) {
      // auto msg = boost::format{"%1% has no atom elements"} %
      // xpath_to_atomArray; throw RDKit::FileParseException{msg.str()};
      throw cml::MandatoryElementNotFoundError{__func__};
    }

    auto conformer = std::make_unique<RDKit::Conformer>(num_atoms);
    unsigned atom_idx = 0u;
    for (const auto &atomitr : *atomArray) {
      if (atomitr.first == "<xmlattr>") {
        for (const auto &i : atomitr.second) {
          BOOST_LOG(rdInfoLog)
              << boost::format{"%1%/@%2% (= \"%3%\") is ignored"} %
                     xpath_to_atomArray % i.first % i.second.data()
              << std::endl;
        }
        continue;
      }

      if (atomitr.first != "atom") {
        BOOST_LOG(rdInfoLog) << boost::format{"%1%/%2% is ignored"} %
                                    xpath_to_atomArray % atomitr.first
                             << std::endl;
        continue;
      }

      // parse_atom(xpath_to_atomArray + "/" + atomitr.first, atomitr.second,
      // atom_idx++);

      const auto &atom_node = atomitr.second;

      // http://www.xml-cml.org/convention/molecular#atom-elementType
      // > An atom MUST have an elementType attribute
      try {
        const auto elementType =
            atom_node.get<std::string>("<xmlattr>.elementType");

        // http://www.xml-cml.org/schema/schema3/schema.xsd
        // > Du. ("dummy") This does not correspond to a "real" atom and can
        // support a > point in space or within a chemical graph. > R.
        // ("R-group") This indicates that an atom or group of atoms could be >
        // attached at this point.
        auto atom = (elementType == "Du" || elementType == "R")
                        ? std::make_unique<RDKit::Atom>()
                        : std::make_unique<RDKit::Atom>(elementType);
      } catch (const boost::property_tree::ptree_bad_path &) {
        auto msg = boost::format{"%1%/@elementType is missing"} % xpath_to_atom;
        throw cml::CMLAttributeError{msg.str()};
      }
    }
  }

  return mol;
}

}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit
