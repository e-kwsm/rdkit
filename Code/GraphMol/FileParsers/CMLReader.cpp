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

#include "FileParsers.h"
#include "CMLReader.h"
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/Invariant.h>

#ifndef __clang__version__    // FIXME
#pragma GCC diagnostic pop    // FIXME
#else                         // FIXME
#pragma clang diagnostic pop  // FIXME
#endif                        // FIXME

// using namespace std::literals::string_literals;  // XXX

namespace RDKit {
namespace v2 {
namespace FileParsers {
CMLSupplier::CMLSupplier(std::unique_ptr<std::istream> &&p_istream,
                         const CMLFileParserParams &params) noexcept(false)
    : p_istream{std::move(p_istream)}, params{params} {
  PRECONDITION(this->p_istream, "bad stream");

  try {
    boost::property_tree::read_xml(*this->p_istream, pt);
  } catch (const boost::property_tree::xml_parser_error &e) {
    throw RDKit::FileParseException{boost::diagnostic_information(e)};
  }

  if (pt.size() > 1u) {
    throw RDKit::FileParseException{"XML MUST NOT have multiple roots"};
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
          << boost::format{"/%1%/molecule is not found"} % root.first
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
    throw FileParseException("EOF hit.");
  }
  auto tmp = parse_molecule_node(molecule_node.value());
  molecule_node = std::nullopt;  // XXX
  return tmp;
}

std::unique_ptr<RWMol> CMLSupplier::parse_molecule_node(
    const boost::property_tree::ptree &node) {
  std::unique_ptr<RWMol> mol;
  return mol;
}

}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit