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

#pragma once

#include <iosfwd>
#include <memory>
#include <string>

// #include <regex>
// #include <unordered_map>
// #include <unordered_set>

// #include <boost/lexical_cast.hpp>
// #include <boost/property_tree/ptree.hpp>
// #include <boost/format.hpp>

#include <RDGeneral/export.h>
// #include <RDGeneral/FileParseException.h>

namespace RDKit {
class RWMol;

namespace v2 {
namespace FileParsers {
struct RDKIT_FILEPARSERS_EXPORT CMLFileParserParams {
  bool sanitize = true;  ///< sanitize the molecule after building it
  bool removeHs = true;  ///< remove Hs after constructing the molecule
};

class RDKIT_FILEPARSERS_EXPORT CMLSupplier {
 public:
  CMLSupplier() = delete;
  CMLSupplier(std::unique_ptr<std::istream> &&p);
  CMLSupplier(const std::string &fileName);
  ~CMLSupplier();
  CMLSupplier(const CMLSupplier &) = delete;
  CMLSupplier &operator=(const CMLSupplier &) = delete;

  void init();
  void reset();
  std::unique_ptr<RWMol> next();
  void close();

 private:
  std::unique_ptr<std::istream> p_istream;
};
}  // namespace FileParsers
}  // namespace v2

#if 0
class Atom;
class Conformer;
class RWMol;

namespace v2 {
namespace FileParsers {
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromCMLDataStream(
    std::istream &inStream,
    const CMLFileParserParams &params = {}) noexcept(false);
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromCMLBlock(
    const std::string &cmlBlock,
    const CMLFileParserParams &params = {}) noexcept(false);
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromCMLFile(
    const std::string &fName,
    const CMLFileParserParams &params = {}) noexcept(false);

// lexical translator that throws exception if a string is not convertible to
// numeric type `T`
template <typename T>
class LexicalTranslator {
 public:
  using internal_type = boost::property_tree::ptree::data_type;
  using external_type = T;
  static_assert(std::is_arithmetic<T>::value);

  LexicalTranslator(std::string path) : path{std::move(path)} {}
  [[nodiscard]] external_type get_value(const internal_type &s) const
      noexcept(false) {
    try {
      return boost::lexical_cast<external_type>(s);
    } catch (const boost::bad_lexical_cast &) {
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
  // @param molecule_xpath XPath for messages
  // @param molecule_node
  CMLMoleculeParser(std::string molecule_xpath,
                    const boost::property_tree::ptree &molecule_node);
  ~CMLMoleculeParser() = default;
  CMLMoleculeParser(const CMLMoleculeParser &) = delete;

  CMLMoleculeParser &operator=(const CMLMoleculeParser &) = delete;

  std::unique_ptr<RWMol> parse(
      const v2::FileParsers::CMLFileParserParams &params);

 private:
  void check_molecule();
  void parse_atomArray(
      const std::string &xpath_to_atomArray,
      const boost::property_tree::ptree &atomArray) noexcept(false);
  void parse_atom(std::string xpath_to_atom,
                  const boost::property_tree::ptree &atom_node,
                  unsigned idx) noexcept(false);
  void parse_bondArray(
      const std::string &xpath_to_bondArray,
      const boost::property_tree::ptree &bondArray) noexcept(false);
  void parse_bond(std::string xpath_to_bond,
                  const boost::property_tree::ptree &bond_node) noexcept(false);

  std::string xpath() const {
    return molecule_xpath +
           (molecule_id ? ("[@id=\"" + *molecule_id + "\"]") : "");
  }
  template <typename S>
  boost::format xpath(const S &path) const {
    return boost::format{"%1%/%2%"} % xpath() % path;
  }

  template <typename T>
  static boost::optional<T> get_attribute_optionally(
      const boost::property_tree::ptree &pt, const std::string &name,
      const std::string &parent) noexcept(false) {
    return pt.get_optional<T>("<xmlattr>." + name,
                              LexicalTranslator<T>{parent + "/@" + name});
  }

  // http://www.xml-cml.org/convention/molecular#molecule-id
  // http://www.xml-cml.org/convention/molecular#atom-id
  // > The value of the id attribute MUST start with a letter,
  // > and MUST only contain letters, numbers, dot, hyphen or underscore.
  static bool is_valid_id(const std::string &id) {
    return std::regex_match(id, std::regex{"^[A-Za-z][A-Za-z0-9._-]*$"});
  }

  void check_hydrogenCount();

 private:
  const std::string molecule_xpath;
  const boost::property_tree::ptree molecule_node;
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
}  // namespace FileParsers
}  // namespace v2
#endif
}  // namespace RDKit
