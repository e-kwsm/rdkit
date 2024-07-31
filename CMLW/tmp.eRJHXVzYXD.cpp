#include "CMLWriter.h"
#include <fstream>

#include <iostream>

#include <boost/format.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <RDGeneral/BadFileException.h>
#include <RDGeneral/FileParseException.h>
#include <GraphMol/RWMol.h>

namespace RDKit {
CMLWriter::CMLWriter(const std::string &fileName)
    : CMLWriter{std::make_unique<std::ofstream>(fileName)} {}

CMLWriter::CMLWriter(std::unique_ptr<std::ostream> &&p)
    : p_ostream{std::move(p)} {
  PRECONDITION(p_ostream, "null stream");
  if (p_ostream->bad()) {
    throw FileParseException("Bad output stream");
  }
  auto &cml = pt.add("cml", "");
  cml.put("<xmlattr>.xmlns", "http://www.xml-cml.org/schema");
  cml.put("<xmlattr>.xmlns:convention", "http://www.xml-cml.org/convention/");
  cml.put("<xmlattr>.convention", "convention:molecular");
}

CMLWriter::~CMLWriter() { write(); }

void CMLWriter::write() const {
  PRECONDITION(p_ostream, "no output stream");
  boost::property_tree::write_xml(
      *p_ostream, pt,
      boost::property_tree::xml_writer_make_settings<std::string>(' ', 2));
}

void CMLWriter::add(const ROMol &mol, int confId) {
  PRECONDITION(p_ostream, "no output stream");
  RWMol rwmol{mol};
  auto &molecule_node = pt.add("cml.molecule", "");
  molecule_node.put("<xmlattr>.id", boost::format{"%1%.%2%"} % "m" % confId);

  constexpr auto atom_id_prefix = "a";

  const Conformer *conf = nullptr;
  if (rwmol.getNumConformers()) {
    conf = &rwmol.getConformer(confId);
  }
  auto &atomArray = molecule_node.put("atomArray", "");
  for (unsigned i = 0u, nAtoms = rwmol.getNumAtoms(); i < nAtoms; i++) {
    auto &atom = atomArray.add("atom", "");
    const auto &a = rwmol.getAtomWithIdx(i);

    atom.put("<xmlattr>.id", boost::format{"%1%%2%"} % atom_id_prefix % i);
    if (a->getAtomicNum()) {
      atom.put("<xmlattr>.elementType", a->getSymbol());
    } else {
      atom.put("<xmlattr>.elementType", "Du");  // dummy
    }

    if (conf != nullptr) {
      const auto &pos = conf->getAtomPos(i);
      boost::format xyz_fmt{"%.6f"};

      if (!conf->is3D()) {
        atom.put("<xmlattr>.x2", xyz_fmt % pos.x);
        atom.put("<xmlattr>.y2", xyz_fmt % pos.y);
      } else {
        atom.put("<xmlattr>.x3", xyz_fmt % pos.x);
        atom.put("<xmlattr>.y3", xyz_fmt % pos.y);
        atom.put("<xmlattr>.z3", xyz_fmt % pos.z);
      }
    }
  }

  auto &bondArray = molecule_node.put("bondArray", "");
  bondArray.add("atom", "");
}
}  // namespace RDKit
