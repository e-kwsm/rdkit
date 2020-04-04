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

#include "./FileParsers.h"  // FIXME

#include <regex>

#include <boost/exception/diagnostic_information.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <GraphMol/FileParsers/FileParsers.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/FileParseException.h>

#ifndef __clang__version__    // FIXME
#pragma GCC diagnostic pop    // FIXME
#else                         // FIXME
#pragma clang diagnostic pop  // FIXME
#endif                        // FIXME

#if 0
#include <iostream>
#include <memory>
#include <unordered_map>

#include <Geometry/point.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/RWMol.h>
#include <RDBoost/python.h>
#include <RDGeneral/export.h>
#endif

// using namespace std::literals::string_literals;

namespace {

#if 0
[[deprecated]] bool ends_with(const std::string& str,
                              const std::string& suffix) {
  if (str.size() < suffix.size()) {
    return false;
  }
  return std::equal(suffix.rbegin(), suffix.rend(), str.rbegin());
}

namespace loglevel {
struct [[deprecated]] Error {};
}  // namespace loglevel

template <typename E>
[[deprecated]] void handler(const std::string& msg) {
  throw E{msg};
}

template <>
[[deprecated]] void handler<loglevel::Error>(const std::string& msg) {
  BOOST_LOG(rdErrorLog) << msg;
}
#endif

#if 0
// @return namespace prefix if `str` is prefixed or not-prefixed `target`
[[deprecated]] boost::optional<std::string> xmlnsXXX(const std::string& target,
                                                     const std::string& str) {
  if (str == target) {  // not prefixed
    return std::string{""};
  }
  std::smatch m;
  if (std::regex_match(
          str, m,
          std::regex{(boost::format{R"((?:(\w*):)?%1%)"} % target).str()})) {
    return {m.str(1)};
  }
  return {};
}

[[deprecated]] std::unordered_map<std::string, std::string> xmlattr(
    const boost::property_tree::ptree& pt) {
  std::unordered_map<std::string, std::string> r;
  for (const auto& it : pt) {
    if (it.first == "<xmlattr>") {
      for (const auto& i : it.second) {
        r[i.first] = i.second.data();
      }
    }
  }
  return r;
}
#endif

}  // namespace

namespace RDKit {
namespace {

// lexical translator that throws exception if a string is not convertible to
// numeric type `T`
template <typename T>
struct LexicalTranslator {
  using internal_type = boost::property_tree::ptree::data_type;
  using external_type = T;
  static_assert(std::is_arithmetic<T>::value, "");

  LexicalTranslator(std::string path) : path{std::move(path)} {}
  external_type get_value(const internal_type& s) const noexcept(false) {
    try {
      return boost::lexical_cast<external_type>(s);
    } catch (const boost::bad_lexical_cast&) {
      auto msg =
          boost::format{"%1% (= \"%2%\") is not convertible to numeric"} %
          path % s;
      throw RDKit::FileParseException{msg.str()};
    }
  }

 private:
  const std::string path;
};

class CMLMoleculeParser {
 public:
  CMLMoleculeParser() = delete;
  // @param molecule_xpath XPath for messages
  // @param molecule_node
  CMLMoleculeParser(std::string molecule_xpath,
                    const boost::property_tree::ptree& molecule_node);
  ~CMLMoleculeParser() = default;

  RDKit::RWMol* parse(bool sanitize, bool removeHs);

 private:
  void parse_atomArray(
      const std::string& xpath_to_atomArray,
      const boost::property_tree::ptree& atomArray) noexcept(false);
  std::pair<std::string, std::unique_ptr<RDKit::Atom>> parse_atom(
      const std::string& xpath_to_atom, const boost::property_tree::ptree& atom,
      unsigned idx) noexcept(false);
  void parse_bondArray(
      const std::string& xpath_to_bondArray,
      const boost::property_tree::ptree& bondArray) noexcept(false);
  std::unique_ptr<RDKit::Bond> parse_bond(
      const std::string& xpath_to_bond,
      const boost::property_tree::ptree& bond);

  std::string xpath() const {
    return molecule_xpath +
           (molecule_id ? ("[@id=\"" + *molecule_id + "\"]") : "");
  }
  template <typename S>
  boost::format xpath(const S& path) const {
    return boost::format{"%1%/%2%"} % xpath() % path;
  }

  template <typename T>
  boost::optional<T> get_attribute_optionally(
      const boost::property_tree::ptree& pt, const std::string& name,
      const std::string& parent) const noexcept(false) {
    return pt.get_optional<T>("<xmlattr>." + name,
                              LexicalTranslator<T>{parent + "/@" + name});
  }

#if 0
  template <typename S>
  boost::format put_missing(const S& path) const {
    return boost::format{"%1% is missing\n"} % path;
  }
#endif

  // http://www.xml-cml.org/convention/molecular#molecule-id
  // http://www.xml-cml.org/convention/molecular#atom-id
  // > The value of the id attribute MUST start with a letter,
  // > and MUST only contain letters, numbers, dot, hyphen or underscore.
  bool is_valid_id(const std::string& id) const {
    return std::regex_match(id, std::regex{"^[A-Za-z][A-Za-z0-9._-]*"});
  }

  void check_hydrogenCount();

 private:
  const boost::property_tree::ptree molecule_node;
  const std::string molecule_xpath;
  const boost::optional<std::string> molecule_id;
  const boost::optional<int> formalCharge;
  const boost::optional<unsigned> spinMultiplicity;

  std::unique_ptr<RDKit::RWMol> molecule;
  std::unique_ptr<RDKit::Conformer> conformer;

  std::unordered_map<std::string, std::unique_ptr<RDKit::Atom>> id_atoms;
  std::unordered_map<std::string, boost::optional<unsigned>> id_hydrogenCounts;
  std::unordered_set<std::string> bond_ids;
};
}  // namespace

CMLMoleculeParser::CMLMoleculeParser(
    std::string molecule_xpath,
    const boost::property_tree::ptree& molecule_node)
    : molecule_node{molecule_node},
      molecule_xpath{std::move(molecule_xpath)},
      molecule_id{molecule_node.get_optional<std::string>("<xmlattr>.id")},
      formalCharge{get_attribute_optionally<int>(molecule_node, "formalCharge",
                                                 this->molecule_xpath)},
      spinMultiplicity{get_attribute_optionally<unsigned>(
          molecule_node, "spinMultiplicity", this->molecule_xpath)} {
  // ignore invalid //molecule/@id
  if (!molecule_id) {
    // http://www.xml-cml.org/convention/molecular#molecule-id
    // > A molecule element MUST have an id attribute
    BOOST_LOG(rdWarningLog)
        << boost::format{"%1%/@id is missing"} % this->molecule_xpath
        << std::endl;
  } else if (!is_valid_id(*molecule_id)) {
    BOOST_LOG(rdWarningLog) << boost::format{"%1%/@id (= \"%2%\") is invalid"} %
                                   this->molecule_xpath % *molecule_id
                            << std::endl;
  }

  if (!formalCharge) {
    // http://www.xml-cml.org/convention/molecular#molecule-formalCharge
    // > A molecule SHOULD have a formalCharge attribute specified
    BOOST_LOG(rdWarningLog)
        << boost::format{"%1% is missing"} % xpath("@formalCharge")
        << std::endl;
  }

  if (!spinMultiplicity) {
    // http://www.xml-cml.org/convention/molecular#molecule-spinMultiplicity
    // > A molecule SHOULD have a spinMultiplicity attribute specified.
    BOOST_LOG(rdWarningLog)
        << boost::format{"%1% is missing"} % xpath("@spinMultiplicity")
        << std::endl;
  } else if (*spinMultiplicity == 0u) {
    auto msg = boost::format{"%1% is zero"} % xpath("@spinMultiplicity");
    throw RDKit::FileParseException{msg.str()};
  }

  molecule = std::make_unique<RDKit::RWMol>();
}

RDKit::RWMol* CMLMoleculeParser::parse(bool sanitize, bool removeHs) {
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
    parse_atomArray(xpath("atomArray").str(), *aa);
  }

  // http://www.xml-cml.org/convention/molecular#molecule-bondArray
  // > A molecule MAY contain a single bondyArray child provided that it does
  // > not contain child molecules.
  if (auto ba = get_array("bondArray")) {
    parse_bondArray(xpath("bondArray").str(), *ba);
  }

  if (sanitize) {
    if (removeHs) {
      MolOps::removeHs(*molecule, false, false);
    } else {
      MolOps::sanitizeMol(*molecule);
    }
    check_hydrogenCount();
  } else {
    // molecule->updatePropertyCache(false);
  }
  // MolOps::assignChiralTypesFrom3D(*molecule);

  return molecule.release();
}

void CMLMoleculeParser::parse_atomArray(
    const std::string& xpath_to_atomArray,
    const boost::property_tree::ptree& atomArray) noexcept(false) {
  const unsigned num_atoms = atomArray.count("atom");
  if (num_atoms == 0u) {
    // http://www.xml-cml.org/convention/molecular#atomArray-element
    // > An atomArray element MUST contain at least one child atom element.
    auto msg = boost::format{"%1% has no atom"} % xpath_to_atomArray;
    throw RDKit::FileParseException{msg.str()};
  }
  conformer = std::make_unique<RDKit::Conformer>(num_atoms);

  unsigned atom_idx = 0u;
  for (const auto& atomitr : atomArray) {
    if (atomitr.first == "<xmlattr>") {
      for (const auto& i : atomitr.second) {
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

    auto id_atom = parse_atom(xpath_to_atomArray + "/" + atomitr.first,
                              atomitr.second, atom_idx++);
    id_atoms[id_atom.first] = std::move(id_atom.second);
  }

  int sum_formalCharge = 0;
  unsigned overall_spinMultiplicity = 1u;
  for (const auto& a : id_atoms) {
    sum_formalCharge += a.second->getFormalCharge();
    overall_spinMultiplicity += a.second->getNumRadicalElectrons();
  }

  if (formalCharge && *formalCharge != sum_formalCharge) {
    auto msg = boost::format{"%1% (= %2%) is not equal to sum of %3% (= %4%)"} %
               xpath("@formalCharge") % *formalCharge %
               "../atomArray/atom/@formalCharge" % sum_formalCharge;
    throw RDKit::FileParseException{msg.str()};
  }

  if (spinMultiplicity && *spinMultiplicity != overall_spinMultiplicity) {
    auto msg =
        boost::format{
            "%1% (= %2%) is not equal to "
            "%3% = %4% (sum of %5%) - %6% (number of atoms)"} %
        xpath("@spinMultiplicity") % *spinMultiplicity %
        overall_spinMultiplicity % (overall_spinMultiplicity + num_atoms) %
        "../atomArray/atom/@spinMultiplicity" % num_atoms;
    throw RDKit::FileParseException{msg.str()};
  }

  for (auto&& p : id_atoms) {
    molecule->addAtom(p.second.release(), true, true);
  }
}

std::pair<std::string, std::unique_ptr<RDKit::Atom>>
CMLMoleculeParser::parse_atom(const std::string& xpath_to_atom,
                              const boost::property_tree::ptree& atom,
                              unsigned idx) noexcept(false) {
  const auto id = atom.get_optional<std::string>("<xmlattr>.id");
  if (!id) {
    // http://www.xml-cml.org/convention/molecular#atom-id
    // > An atom MUST have an id attribute
    auto msg = boost::format{"%1%/@id is missing"} % xpath_to_atom;
    throw RDKit::FileParseException{msg.str()};
  } else if (!is_valid_id(*id)) {
    auto msg =
        boost::format{"%1%/@id (= \"%2%\") is invalid"} % xpath_to_atom % *id;
    throw RDKit::FileParseException{msg.str()};
  } else if (id_atoms.count(*id)) {
    // http://www.xml-cml.org/convention/molecular#atom-id
    // > The value of the id MUST be unique amongst the atoms within the eldest
    // > containing molecule.
    auto msg = boost::format{"%1%/@id (= \"%2%\") is not unique"} %
               xpath_to_atom % *id;
    throw RDKit::FileParseException{msg.str()};
  }
  const auto xpath_to_atom_w_id =
      (boost::format{"%1%[@id=\"%2%\"]"} % xpath_to_atom % *id).str();

  const auto elementType =
      atom.get_optional<std::string>("<xmlattr>.elementType");
  if (!elementType) {
    // http://www.xml-cml.org/convention/molecular#atom-elementType
    // > An atom MUST have an elementType attribute
    auto msg =
        boost::format{"%1%/@elementType is missing"} % xpath_to_atom_w_id;
    throw RDKit::FileParseException{msg.str()};
  }

#if 0
  BOOST_LOG(rdDebugLog) << xpath_to_atom_w_id << elementType << '\n';
#endif

  // http://www.xml-cml.org/schema/schema3/schema.xsd
  // > Du. ("dummy") This does not correspond to a "real" atom and can support a
  // > point in space or within a chemical graph.
  // > R. ("R-group") This indicates that an atom or group of atoms could be
  // > attached at this point.
  auto a = (*elementType == "Du" || *elementType == "R")
               ? std::make_unique<RDKit::Atom>()
               : std::make_unique<RDKit::Atom>(*elementType);

  // http://www.xml-cml.org/convention/molecular#atom-isotopeNumber
  // > An atom MAY have an isotopeNumber attribute.
  if (auto isotopeNumber = atom.get_optional<unsigned>(
          "<xmlattr>.isotopeNumber",
          LexicalTranslator<unsigned>{xpath_to_atom_w_id +
                                      "/@isotopeNumber"})) {
    a->setIsotope(*isotopeNumber);
  }

  // http://www.xml-cml.org/convention/molecular#atom-formalCharge
  // > An atom MAY have a formalCharge attribute.
  if (auto formalCharge = atom.get_optional<int>(
          "<xmlattr>.formalCharge",
          LexicalTranslator<int>{xpath_to_atom_w_id + "/@formalCharge"})) {
    a->setFormalCharge(*formalCharge);
  }

  // http://www.xml-cml.org/convention/molecular#atom-spinMultiplicity
  // > An atom MAY have a spinMultiplicity attribute.
  if (auto spinMultiplicity = atom.get_optional<unsigned>(
          "<xmlattr>.spinMultiplicity",
          LexicalTranslator<unsigned>{xpath_to_atom_w_id +
                                      "/@spinMultiplicity"})) {
    if (*spinMultiplicity == 0u) {
      auto msg =
          boost::format{"%1%/@spinMultiplicity is zero"} % xpath_to_atom_w_id;
      throw RDKit::FileParseException{msg.str()};
    }
    a->setNumRadicalElectrons(*spinMultiplicity - 1u);
  }

  const auto hydrogenCount = atom.get_optional<unsigned>(
      "<xmlattr>.hydrogenCount",
      LexicalTranslator<unsigned>{xpath_to_atom_w_id + "/@hydrogenCount"});
  id_hydrogenCounts[xpath_to_atom_w_id] = hydrogenCount;

  // http://www.xml-cml.org/convention/molecular#atom-x3
  // > An atom MAY have an x3 attribute, the value of which is the x coordinate
  // > of a 3 dimensional object. The units are Angstrom and the axis system is
  // > always right handed.
  const auto x3 = atom.get_optional<double>(
      "<xmlattr>.x3", LexicalTranslator<double>{xpath_to_atom_w_id + "/@x3"});
  const auto y3 = atom.get_optional<double>(
      "<xmlattr>.y3", LexicalTranslator<double>{xpath_to_atom_w_id + "/@y3"});
  const auto z3 = atom.get_optional<double>(
      "<xmlattr>.z3", LexicalTranslator<double>{xpath_to_atom_w_id + "/@z3"});

  // http://www.xml-cml.org/convention/molecular#atom-x2
  // > An atom MAY have an x2 attribute, the value of which is used for
  // > displaying the object in 2 dimensions. This is unrelated to the 3-D
  // > coordinates for the object.
  const auto x2 = atom.get_optional<double>(
      "<xmlattr>.x2", LexicalTranslator<double>{xpath_to_atom_w_id + "/@x2"});
  const auto y2 = atom.get_optional<double>(
      "<xmlattr>.y2", LexicalTranslator<double>{xpath_to_atom_w_id + "/@x2"});

  if (x3 || y3 || z3) {
    if (x2 || y2) {
      auto msg =
          boost::format{
              "%1% has both of 3D (any of @x3, @y3 and @z3) "
              "and 2D (any of @x2 and @y2) coordinate attributes"} %
          xpath_to_atom_w_id;
      throw RDKit::FileParseException{msg.str()};
    }

    if (!(x3 && y3 && z3)) {
      // http://www.xml-cml.org/convention/molecular#atom-x3
      // > If a x3 attribute is present there MUST also be a y3 and z3 present.
      auto msg =
          boost::format{"%1% does not have all of x3, y3 and z3 attributes"} %
          xpath_to_atom_w_id;
      throw RDKit::FileParseException{msg.str()};
    }
    RDGeom::Point3D r{*x3, *y3, *z3};
    BOOST_LOG(rdDebugLog) << xpath_to_atom_w_id << ' ' << r << std::endl;
    conformer->setAtomPos(idx, r);
  } else if (x2 || y2) {
    if (!(x2 && y2)) {
      // http://www.xml-cml.org/convention/molecular#atom-x2
      // > If a x2 attribute is present there MUST also be a y2 attribute.
      auto msg =
          boost::format{"%1% does not have both of x2 and y2 attributes"} %
          xpath_to_atom_w_id;
      throw RDKit::FileParseException{msg.str()};
    }
    RDGeom::Point3D r{*x2, *y2, 0.0};
    BOOST_LOG(rdDebugLog) << xpath_to_atom_w_id << ' ' << r << std::endl;
    conformer->setAtomPos(idx, r);
  } else {
    BOOST_LOG(rdInfoLog)
        << boost::format{"%1% does not have geometrical info"} %
               xpath_to_atom_w_id
        << std::endl;
  }

  return {*id, std::move(a)};
}

void CMLMoleculeParser::parse_bondArray(
    const std::string& xpath_to_bondArray,
    const boost::property_tree::ptree& bondArray) noexcept(false) {
  if (bondArray.count("bond") == 0u) {
    // http://www.xml-cml.org/convention/molecular#bondArray-element
    // > A bondArray element MUST contain at least one child bond element.
    auto msg = boost::format{"%1% has no bond"} % xpath_to_bondArray;
    throw RDKit::FileParseException{msg.str()};
  }

  for (const auto& bonditr : bondArray) {
    if (bonditr.first == "<xmlattr>") {
#if 0
      for (const auto& i : bonditr.second) {
        std::cerr << __FILE__ << ':' << __LINE__ << '\t';
        std::cerr << boost::format{"%1%/@%2% = \"%3%\""} % xpath_to_bondArray %
                         i.first % i.second.data();
        std::cerr << '\n';
      }
#endif
      continue;
    }

    if (bonditr.first != "bond") {
      BOOST_LOG(rdInfoLog) << boost::format{"%1%/%2% is ignored"} %
                                  xpath_to_bondArray % bonditr.first
                           << std::endl;
      continue;
    }

    if (auto b = parse_bond(xpath_to_bondArray + "/" + bonditr.first,
                            bonditr.second)) {
      molecule->addBond(b.release(), true);
    }
  }
}

std::unique_ptr<RDKit::Bond> CMLMoleculeParser::parse_bond(
    const std::string& xpath_to_bond, const boost::property_tree::ptree& bond) {
  // http://www.xml-cml.org/convention/molecular#bond-id
  // > It is RECOMMENDED that a bond has an id attribute so that it can be
  // > referenced.
  const auto id = bond.get_optional<std::string>("<xmlattr>.id");
  if (!id) {
    BOOST_LOG(rdInfoLog) << boost::format{"%1%/@id is missing"} % xpath_to_bond
                         << std::endl;
  } else if (!is_valid_id(*id)) {
    BOOST_LOG(rdWarningLog)
        << boost::format{"%1%/@id (= \"%2%\") is invalid"} % xpath_to_bond % id
        << std::endl;
  } else if (!bond_ids.insert(*id).second) {
    // http://www.xml-cml.org/convention/molecular#bond-id
    // > The id of a bond MUST be unique amongst the bonds of the eldest
    // > containing molecule.
    auto msg = boost::format{"%1%/@id (= \"%2%\") is not unique"} %
               xpath_to_bond % *id;
    throw RDKit::FileParseException{msg.str()};
  }
  const auto xpath_to_bond_w_id =
      xpath_to_bond + (id ? ("[@id=\"" + *id + "\"]") : "");

  const auto order = bond.get_optional<std::string>("<xmlattr>.order");
  if (!order) {
    // http://www.xml-cml.org/convention/molecular#bond-order
    // > A bond MUST have an order attribute
    auto msg = boost::format{"%1%/@order is missing"} % xpath_to_bond_w_id;
    throw RDKit::FileParseException{msg.str()};
  }

  RDKit::Bond::BondType bt = RDKit::Bond::UNSPECIFIED;
  if (*order == "1" || *order == "S") {
    bt = RDKit::Bond::SINGLE;
  } else if (*order == "2" || *order == "D") {
    bt = RDKit::Bond::DOUBLE;
  } else if (*order == "3" || *order == "T") {
    bt = RDKit::Bond::TRIPLE;
  } else if (*order == "A") {
    bt = RDKit::Bond::AROMATIC;
  } /* RDKit extension */ else if (*order == "4") {
    bt = RDKit::Bond::QUADRUPLE;
  } else if (*order == "5") {
    bt = RDKit::Bond::QUINTUPLE;
  } else if (*order == "6") {
    bt = RDKit::Bond::HEXTUPLE;
  } else if (*order == "1.5") {
    bt = RDKit::Bond::ONEANDAHALF;
  } else if (*order == "2.5") {
    bt = RDKit::Bond::TWOANDAHALF;
  } else if (*order == "3.5") {
    bt = RDKit::Bond::TWOANDAHALF;
  } else if (*order == "4.5") {
    bt = RDKit::Bond::FOURANDAHALF;
  } else if (*order == "5.5") {
    bt = RDKit::Bond::FIVEANDAHALF;
  } else {
    BOOST_LOG(rdWarningLog)
        << boost::format{"%1%/@order (= \"%2%\") is unrecognizable"} %
               xpath_to_bond_w_id % *order
        << std::endl;
  }

  // http://www.xml-cml.org/convention/molecular#bond-atomRefs2
  // > A bond MUST have a atomRefs2 attribute, the value of which MUST be the
  // > space separated ids of two different atoms which MUST be in the same
  // > molecule.
  const auto atomRefs2 = bond.get_optional<std::string>("<xmlattr>.atomRefs2");
  if (!atomRefs2) {
    auto msg = boost::format{"%1%/@atomRefs2 is missing"} % xpath_to_bond_w_id;
    throw RDKit::FileParseException{msg.str()};
  }
  std::istringstream iss{*atomRefs2};
  std::string id_bgn, id_end, extra;
  if (!(iss >> id_bgn >> id_end)) {
    auto msg =
        boost::format{
            "%1%/@atomRefs2 (= \"%2%\") does not have two ids "
            "separated by space"} %
        xpath_to_bond_w_id % *atomRefs2;
    throw RDKit::FileParseException{msg.str()};
  } else if (iss >> extra) {
    auto msg =
        boost::format{"%1%/@atomRefs2 (= \"%2%\") has three or more ids"} %
        xpath_to_bond_w_id % *atomRefs2;
    throw RDKit::FileParseException{msg.str()};
  }

  if (id_bgn == id_end) {
    auto msg = boost::format{"%1%/@atomRefs2 (= \"%2%\") is self-bond"} %
               xpath_to_bond_w_id % *atomRefs2;
    throw RDKit::FileParseException{msg.str()};
  }

  auto distance = [&](const auto& atom_id) {
    auto i = id_atoms.find(atom_id);
    if (i == id_atoms.end()) {
      auto msg =
          boost::format{
              "%1%/@atomRefs2 (= \"%2%\") refers to non-existing "
              "../../../atomArray/atom[@id=\"%3%\"]"} %
          xpath_to_bond_w_id % *atomRefs2 % atom_id;
      throw RDKit::FileParseException{msg.str()};
    }
    return std::distance(id_atoms.begin(), i);
  };

  const auto index_bgn = distance(id_bgn);
  const auto index_end = distance(id_end);

  auto b = std::make_unique<RDKit::Bond>(bt);
  b->setBeginAtomIdx(index_bgn);
  b->setEndAtomIdx(index_end);

  // TODO
  const auto bondStereo = bond.get_optional<std::string>("bondStereo");
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
                 xpath_to_bond_w_id % *bondStereo
          << std::endl;
    }
  }

  b->setOwningMol(*molecule);

  return b;
}

void CMLMoleculeParser::check_hydrogenCount() {
  unsigned i = 0u;
  for (const auto& id_hc : id_hydrogenCounts) {
    const auto& a = (*molecule)[i++];
    a->calcExplicitValence(false);
    unsigned nh = 0u;
    for (const auto neighbor :
         boost::make_iterator_range(molecule->getAtomNeighbors(a))) {
      if ((*molecule)[neighbor]->getAtomicNum() == 1) {
        nh++;
      }
    }
    BOOST_LOG(rdDebugLog) << *a << " " << nh << "\n";

    const auto& hydrogenCount = id_hc.second;
    if (!hydrogenCount) {
      continue;
    }
    if (*hydrogenCount < nh) {
      BOOST_LOG(rdWarningLog) <<
          boost::format{
              "%1%/@hydrogenCount (= %2%) is less than "
              "the number of explicitly connected hydrogens (= %3%)"} %
              id_hc.first % *hydrogenCount % nh
          << std::endl;
    } else {
      a->setNumExplicitHs(nh);
    }
  }
}

RDKit::RWMol* PTreeToMol(const boost::property_tree::ptree& pt, bool sanitize,
                         bool removeHs) noexcept(false) {
  if (pt.size() > 1u) {
    throw RDKit::FileParseException{"XML MUST NOT have multiple roots"};
  }
  const auto& root = pt.front();

#if 0
  for (const auto& it : root.second) {
    if (it.first == "<xmlattr>") {  // attribute
      for (const auto& i : it.second) {
        std::cerr << __FILE__ << ':' << __LINE__ << '\t';
        std::cerr << boost::format{"/%1%/@%2% = \"%3%\""} % root.first %
                         i.first % i.second.data();
        std::cerr << '\n';
      }
      if (auto p = xmlnsXXX("molecule", it.first)) {
        std::cerr << __FILE__ << ':' << __LINE__ << '\t';
        std::cerr << *p << '\n';
      }
    }
    if (auto p = xmlnsXXX("molecule", it.first)) {
      std::cerr << __FILE__ << ':' << __LINE__ << '\t';
      std::cerr << *p << '\n';
    }
  }
#endif

#if 0
  if (auto xmlns = root.second.get_optional<std::string>("<xmlattr>.xmlns")) {
    std::cerr << __FILE__ << ':' << __LINE__ << '\t' << *xmlns << '\n';
  }
  if (auto xmlns_cml =
          root.second.get_optional<std::string>("<xmlattr>.xmlns:cml")) {
    std::cerr << __FILE__ << ':' << __LINE__ << '\t' << *xmlns_cml << '\n';
  }
  if (auto xmlns_convention =
          root.second.get_optional<std::string>("<xmlattr>.xmlns:convention")) {
    std::cerr << __FILE__ << ':' << __LINE__ << '\t' << *xmlns_convention
              << '\n';
  }
#endif

  std::unique_ptr<CMLMoleculeParser> p;
  if (root.first == "molecule") {
    p = std::make_unique<CMLMoleculeParser>("/" + root.first, root.second);
  } else {
    // http://www.xml-cml.org/convention/molecular#applying
    // > If the molecular convention is specified on a cml element then that
    // > element MUST have at least one child molecule element that either has
    // > no convention specified or specifies the molecular convention.
    if (root.first != "cml") {
      BOOST_LOG(rdWarningLog)
          << boost::format{"XML root element is %1%"} % root.first << std::endl;
    }
    auto m = root.second.find("molecule");
    // XXX do not dig into 3rd or deeper
    if (m == root.second.not_found()) {
      BOOST_LOG(rdWarningLog)
          << boost::format{"/%1%/molecule is not found"} % root.first
          << std::endl;
      return nullptr;
    }
    // XXX `cml` may have multiple `molecule` elements
    // XXX nested `molecule` elements are not supported
    if (root.second.count("molecule") > 1u) {
      BOOST_LOG(rdWarningLog)
          << boost::format{"/%1% has multiple molecule elements"} % root.first
          << std::endl;
    }
    p = std::make_unique<CMLMoleculeParser>("/" + root.first + "/" + m->first,
                                            m->second);
  }
  return p->parse(sanitize, removeHs);
}

RDKit::RWMol* CMLBlockToMol(const std::string& block, bool sanitize,
                            bool removeHs) noexcept(false) {
  std::istringstream iss{block};
  boost::property_tree::ptree pt;
  try {
    boost::property_tree::read_xml(iss, pt);
  } catch (const boost::property_tree::xml_parser_error& e) {
    throw FileParseException{boost::diagnostic_information(e)};
  }
  return PTreeToMol(pt, sanitize, removeHs);
}

RDKit::RWMol* CMLFileToMol(const std::string& filename, bool sanitize,
                           bool removeHs) noexcept(false) {
  {  // distinguish IOError from other errors, e.g., malformed XML
    std::ifstream ifs{filename};
    if (!ifs || ifs.bad()) {
      throw BadFileException{
          (boost::format{"Bad input file %1%"} % filename).str()};
    }
  }

  boost::property_tree::ptree pt;
  try {
    boost::property_tree::read_xml(filename, pt);
  } catch (const boost::property_tree::xml_parser_error& e) {
    throw FileParseException{boost::diagnostic_information(e)};
  }
  return PTreeToMol(pt, sanitize, removeHs);
}
}  // namespace RDKit
