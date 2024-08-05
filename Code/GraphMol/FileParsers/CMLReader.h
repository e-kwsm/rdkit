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
#include <optional>
#include <string>

#include <boost/property_tree/ptree.hpp>

#include <RDGeneral/FileParseException.h>
#include <RDGeneral/export.h>

namespace RDKit {
class RWMol;

namespace v2 {
namespace FileParsers {
struct RDKIT_FILEPARSERS_EXPORT CMLFileParserParams {
  bool sanitize = true;  ///< sanitize the molecule after building it
  bool removeHs = true;  ///< remove Hs after constructing the molecule
};

namespace cml {
class CMLError : public RDKit::FileParseException {
 public:
  CMLError(const std::string &what) : RDKit::FileParseException{what} {}
};

class XMLMalformedError : public CMLError {
 public:
  XMLMalformedError(const std::string &what) : CMLError{what} {}
};

class MandatoryElementNotFound : public CMLError {
 public:
  MandatoryElementNotFound(const std::string &what) : CMLError{what} {}
};
}  // namespace cml

class RDKIT_FILEPARSERS_EXPORT CMLSupplier {
 public:
  CMLSupplier() = delete;
  CMLSupplier(std::unique_ptr<std::istream> &&p_istream,
              const CMLFileParserParams &params = {}) noexcept(false);
  CMLSupplier(const std::string &fileName,
              const CMLFileParserParams &params = {}) noexcept(false);
  ~CMLSupplier();
  CMLSupplier(const CMLSupplier &) = delete;
  CMLSupplier &operator=(const CMLSupplier &) = delete;

  void reset();
  std::unique_ptr<RWMol> next();
  void close();

 private:
  std::unique_ptr<boost::property_tree::ptree> get_array(
      const boost::property_tree::ptree &node, const std::string &name) const
      noexcept(false);
  std::unique_ptr<RWMol> parse_molecule_node(
      const boost::property_tree::ptree &node);
  std::unique_ptr<RWMol> parse_atom_node(
      const boost::property_tree::ptree &node);
  std::unique_ptr<RWMol> parse_bond_node(
      const boost::property_tree::ptree &node);

  std::unique_ptr<std::istream> p_istream;
  const CMLFileParserParams params;
  boost::property_tree::ptree pt;
  std::optional<boost::property_tree::ptree> molecule_node;
};
}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit
