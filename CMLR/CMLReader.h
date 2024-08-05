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

#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include <boost/property_tree/ptree.hpp>

#include <RDGeneral/FileParseException.h>

namespace RDKit {
class Atom;
class Conformer;
class RWMol;

namespace v2 {
namespace FileParsers {
namespace cml {
class CMLError : public RDKit::FileParseException {
 public:
  CMLError(const std::string &what) : RDKit::FileParseException{what} {}
  ~CMLError() override;
};

class CMLMolecule {
 public:
  CMLMolecule(const boost::property_tree::ptree &molecule_node);
  ~CMLMolecule();

  std::unique_ptr<RWMol> parse();

 private:
  auto get_array(const std::string &name) const;
  void parse_atomArray(const boost::property_tree::ptree &node);
  void parse_atom(const boost::property_tree::ptree &node, unsigned idx);
  void parse_bondArray(const boost::property_tree::ptree &node);
  void parse_bond(const boost::property_tree::ptree &node);

  boost::property_tree::ptree molecule_node;

  std::unique_ptr<RDKit::RWMol> molecule;
  std::unique_ptr<RDKit::Conformer> conformer;
  std::unordered_map<std::string, std::unique_ptr<RDKit::Atom>> id_atom;
  std::unordered_map<std::string, boost::optional<unsigned>> id_hydrogenCount;
  std::unordered_set<std::string> bond_ids;
  int sum_formalCharge = 0;
};
}  // namespace cml
}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit
