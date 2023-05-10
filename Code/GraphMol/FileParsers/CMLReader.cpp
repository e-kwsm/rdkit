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

// using namespace std::literals::string_literals;  // XXX

namespace RDKit {
namespace {

// lexical translator that throws exception if a string is not convertible to
// numeric type `T`
template <typename T>
class LexicalTranslator {
 public:
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
                    const boost::property_tree::ptree& molecule_node)
      : molecule_node{molecule_node},
        molecule_xpath{std::move(molecule_xpath)},
        molecule_id{molecule_node.get_optional<std::string>("<xmlattr>.id")},
        formalCharge{get_attribute_optionally<int>(molecule_node,
                                                   "formalCharge", xpath())},
        spinMultiplicity{get_attribute_optionally<unsigned>(
            molecule_node, "spinMultiplicity", xpath())} {}
  ~CMLMoleculeParser() = default;

  RDKit::RWMol* parse(bool sanitize, bool removeHs);

 private:
  void check_molecule();
  void parse_atomArray(
      const std::string& xpath_to_atomArray,
      const boost::property_tree::ptree& atomArray) noexcept(false);
  void parse_atom(std::string xpath_to_atom,
                  const boost::property_tree::ptree& atom_node,
                  unsigned idx) noexcept(false);
  void parse_bondArray(
      const std::string& xpath_to_bondArray,
      const boost::property_tree::ptree& bondArray) noexcept(false);
  void parse_bond(std::string xpath_to_bond,
                  const boost::property_tree::ptree& bond_node) noexcept(false);

  std::string xpath() const {
    return molecule_xpath +
           (molecule_id ? ("[@id=\"" + *molecule_id + "\"]") : "");
  }
  template <typename S>
  boost::format xpath(const S& path) const {
    return boost::format{"%1%/%2%"} % xpath() % path;
  }

  template <typename T>
  static boost::optional<T> get_attribute_optionally(
      const boost::property_tree::ptree& pt, const std::string& name,
      const std::string& parent) noexcept(false) {
    return pt.get_optional<T>("<xmlattr>." + name,
                              LexicalTranslator<T>{parent + "/@" + name});
  }

  // http://www.xml-cml.org/convention/molecular#molecule-id
  // http://www.xml-cml.org/convention/molecular#atom-id
  // > The value of the id attribute MUST start with a letter,
  // > and MUST only contain letters, numbers, dot, hyphen or underscore.
  static bool is_valid_id(const std::string& id) {
    return std::regex_match(id, std::regex{"^[A-Za-z][A-Za-z0-9._-]*$"});
  }

  void check_hydrogenCount();

 private:
  const boost::property_tree::ptree molecule_node;
  const std::string molecule_xpath;
  const boost::optional<std::string> molecule_id;
  const boost::optional<int> formalCharge;
  const boost::optional<unsigned> spinMultiplicity;

  std::unique_ptr<RDKit::RWMol> molecule = std::make_unique<RDKit::RWMol>();
  std::unique_ptr<RDKit::Conformer> conformer =
      std::make_unique<RDKit::Conformer>();

  std::unordered_map<std::string, std::unique_ptr<RDKit::Atom>> id_atom;
  std::unordered_map<std::string, boost::optional<unsigned>> id_hydrogenCount;
  std::unordered_set<std::string> bond_ids;
  int sum_formalCharge = 0;
};
}  // namespace

void CMLMoleculeParser::check_molecule() {
  // http://www.xml-cml.org/convention/molecular#molecule-id
  // > A molecule element MUST have an id attribute
  // XXX ignore non-existing or invalid //molecule/@id
  if (!molecule_id) {
    BOOST_LOG(rdWarningLog)
        << boost::format{"%1%/@id is missing"} % this->molecule_xpath
        << std::endl;
  } else if (!is_valid_id(*molecule_id)) {
    BOOST_LOG(rdWarningLog) << boost::format{"%1%/@id (= \"%2%\") is invalid"} %
                                   this->molecule_xpath % *molecule_id
                            << std::endl;
  }

  // http://www.xml-cml.org/convention/molecular#molecule-formalCharge
  // > A molecule SHOULD have a formalCharge attribute specified
  if (!formalCharge) {
    BOOST_LOG(rdWarningLog)
        << boost::format{"%1% is missing"} % xpath("@formalCharge")
        << std::endl;
  }

  // http://www.xml-cml.org/convention/molecular#molecule-spinMultiplicity
  // > A molecule SHOULD have a spinMultiplicity attribute specified.
  if (!spinMultiplicity) {
    BOOST_LOG(rdWarningLog)
        << boost::format{"%1% is missing"} % xpath("@spinMultiplicity")
        << std::endl;
  } else if (*spinMultiplicity == 0u) {
    auto msg = boost::format{"%1% is zero"} % xpath("@spinMultiplicity");
    throw RDKit::FileParseException{msg.str()};
  }
}

RDKit::RWMol* CMLMoleculeParser::parse(bool sanitize, bool removeHs) {
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
    parse_atomArray(xpath("atomArray").str(), *aa);
  }

  // http://www.xml-cml.org/convention/molecular#molecule-bondArray
  // > A molecule MAY contain a single bondyArray child provided that it does
  // > not contain child molecules.
  if (auto ba = get_array("bondArray")) {
    parse_bondArray(xpath("bondArray").str(), *ba);
  }

  if (conformer) {
    molecule->addConformer(conformer.release());
  }

  check_hydrogenCount();
  if (sanitize) {
    if (removeHs) {
      MolOps::removeHs(*molecule, false, false);
    } else {
      MolOps::sanitizeMol(*molecule);
    }
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
  // http://www.xml-cml.org/convention/molecular#atomArray-element
  // > An atomArray element MUST contain at least one child atom element.
  if (num_atoms == 0u) {
    auto msg = boost::format{"%1% has no atom elements"} % xpath_to_atomArray;
    throw RDKit::FileParseException{msg.str()};
  }
  conformer->reserve(num_atoms);

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

    parse_atom(xpath_to_atomArray + "/" + atomitr.first, atomitr.second,
               atom_idx++);
  }

  if (formalCharge && *formalCharge != sum_formalCharge) {
    auto msg = boost::format{"%1% (= %2%) is not equal to sum of %3% (= %4%)"} %
               xpath("@formalCharge") % *formalCharge %
               "../atomArray/atom/@formalCharge" % sum_formalCharge;
    throw RDKit::FileParseException{msg.str()};
  }

  // TODO check spinMultiplicity

  for (auto&& p : id_atom) {
    molecule->addAtom(p.second.release(), true, true);
  }
}

void CMLMoleculeParser::parse_atom(std::string xpath_to_atom,
                                   const boost::property_tree::ptree& atom_node,
                                   unsigned idx) noexcept(false) {
  const auto id = atom_node.get_optional<std::string>("<xmlattr>.id");
  if (!id) {
    // http://www.xml-cml.org/convention/molecular#atom-id
    // > An atom MUST have an id attribute
    auto msg = boost::format{"%1%/@id is missing"} % xpath_to_atom;
    throw RDKit::FileParseException{msg.str()};
  } else if (!is_valid_id(*id)) {
    auto msg =
        boost::format{"%1%/@id (= \"%2%\") is invalid"} % xpath_to_atom % *id;
    throw RDKit::FileParseException{msg.str()};
  } else if (id_atom.count(*id)) {
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
    // Ignore x2 and y2 since they are used for displaying the object in 2
    // dimensions
    // http://www.xml-cml.org/convention/molecular#atom-x2
    BOOST_LOG(rdInfoLog)
        << boost::format{"%1% does not have geometrical info"} % xpath_to_atom
        << std::endl;
  }

  id_atom[*id] = std::move(atom);
}

void CMLMoleculeParser::parse_bondArray(
    const std::string& xpath_to_bondArray,
    const boost::property_tree::ptree& bondArray) noexcept(false) {
  if (bondArray.count("bond") == 0u) {
    // http://www.xml-cml.org/convention/molecular#bondArray-element
    // > A bondArray element MUST contain at least one child bond element.
    auto msg = boost::format{"%1% has no bond elements"} % xpath_to_bondArray;
    throw RDKit::FileParseException{msg.str()};
  }

  for (const auto& bonditr : bondArray) {
    if (bonditr.first == "<xmlattr>") {
      continue;
    }

    if (bonditr.first != "bond") {
      BOOST_LOG(rdInfoLog) << boost::format{"%1%/%2% is ignored"} %
                                  xpath_to_bondArray % bonditr.first
                           << std::endl;
      continue;
    }

    parse_bond(xpath_to_bondArray + "/" + bonditr.first, bonditr.second);
  }
}

void CMLMoleculeParser::parse_bond(
    std::string xpath_to_bond,
    const boost::property_tree::ptree& bond_node) noexcept(false) {
  // http://www.xml-cml.org/convention/molecular#bond-id
  // > It is RECOMMENDED that a bond has an id attribute so that it can be
  // > referenced.
  const auto id = bond_node.get_optional<std::string>("<xmlattr>.id");
  if (!id) {
    BOOST_LOG(rdInfoLog) << boost::format{"%1%/@id is missing"} % xpath_to_bond
                         << std::endl;
  } else if (!is_valid_id(*id)) {
    auto msg =
        boost::format{"%1%/@id (= \"%2%\") is invalid"} % xpath_to_bond % id;
    throw RDKit::FileParseException{msg.str()};
  } else if (!bond_ids.insert(*id).second) {
    // http://www.xml-cml.org/convention/molecular#bond-id
    // > The id of a bond MUST be unique amongst the bonds of the eldest
    // > containing molecule.
    auto msg = boost::format{"%1%/@id (= \"%2%\") is not unique"} %
               xpath_to_bond % *id;
    throw RDKit::FileParseException{msg.str()};
  } else {
    xpath_to_bond += "[@id=\"" + *id + "\"]";
  }

  const auto order = bond_node.get_optional<std::string>("<xmlattr>.order");
  if (!order) {
    // http://www.xml-cml.org/convention/molecular#bond-order
    // > A bond MUST have an order attribute
    auto msg = boost::format{"%1%/@order is missing"} % xpath_to_bond;
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
    bt = RDKit::Bond::THREEANDAHALF;
  } else if (*order == "4.5") {
    bt = RDKit::Bond::FOURANDAHALF;
  } else if (*order == "5.5") {
    bt = RDKit::Bond::FIVEANDAHALF;
  } else {
    BOOST_LOG(rdWarningLog)
        << boost::format{"%1%/@order (= \"%2%\") is unrecognizable"} %
               xpath_to_bond % *order
        << std::endl;
  }

  // http://www.xml-cml.org/convention/molecular#bond-atomRefs2
  // > A bond MUST have a atomRefs2 attribute, the value of which MUST be the
  // > space separated ids of two different atoms which MUST be in the same
  // > molecule.
  const auto atomRefs2 =
      bond_node.get_optional<std::string>("<xmlattr>.atomRefs2");
  if (!atomRefs2) {
    auto msg = boost::format{"%1%/@atomRefs2 is missing"} % xpath_to_bond;
    throw RDKit::FileParseException{msg.str()};
  }
  std::istringstream iss{*atomRefs2};
  std::string id_bgn, id_end, extra;
  if (!(iss >> id_bgn >> id_end)) {
    auto msg =
        boost::format{
            "%1%/@atomRefs2 (= \"%2%\") does not have two ids "
            "separated by space"} %
        xpath_to_bond % *atomRefs2;
    throw RDKit::FileParseException{msg.str()};
  } else if (iss >> extra) {
    auto msg =
        boost::format{"%1%/@atomRefs2 (= \"%2%\") has three or more ids"} %
        xpath_to_bond % *atomRefs2;
    throw RDKit::FileParseException{msg.str()};
  }

  if (id_bgn == id_end) {
    auto msg = boost::format{"%1%/@atomRefs2 (= \"%2%\") is self-bond"} %
               xpath_to_bond % *atomRefs2;
    throw RDKit::FileParseException{msg.str()};
  }

  auto distance = [&](const auto& atom_id) {
    auto i = id_atom.find(atom_id);
    if (i == id_atom.end()) {
      auto msg =
          boost::format{
              "%1%/@atomRefs2 (= \"%2%\") refers to non-existing "
              "../../../atomArray/atom[@id=\"%3%\"]"} %
          xpath_to_bond % *atomRefs2 % atom_id;
      throw RDKit::FileParseException{msg.str()};
    }
    return std::distance(id_atom.begin(), i);
  };

  const auto index_bgn = distance(id_bgn);
  const auto index_end = distance(id_end);

  auto bond = std::make_unique<RDKit::Bond>(bt);
  bond->setBeginAtomIdx(index_bgn);
  bond->setEndAtomIdx(index_end);

  // TODO
  const auto bondStereo = bond_node.get_optional<std::string>("bondStereo");
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

void CMLMoleculeParser::check_hydrogenCount() {
  unsigned i = 0u;
  for (const auto& id_hc : id_hydrogenCount) {
    const auto& atom = (*molecule)[i++];
    const auto& hydrogenCount = id_hc.second;
    if (!hydrogenCount) {
      continue;
    }

    atom->calcExplicitValence(false);
    unsigned nh = 0u;
    for (const auto neighbor :
         boost::make_iterator_range(molecule->getAtomNeighbors(atom))) {
      if ((*molecule)[neighbor]->getAtomicNum() == 1) {
        nh++;
      }
    }

    if (*hydrogenCount < nh) {
      BOOST_LOG(rdWarningLog) <<
          boost::format{
              "%1%/@hydrogenCount (= %2%) is less than "
              "the number of explicitly connected hydrogens (= %3%)"} %
              id_hc.first % *hydrogenCount % nh
          << std::endl;
    } else {
      atom->setNumExplicitHs(nh);
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
