#include "CMLWriter.h"
#include <fstream>

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
  auto &cml = tree.add("cml", "");
  cml.put("<xmlattr>.xmlns", "http://www.xml-cml.org/schema");
  cml.put("<xmlattr>.xmlns:convention", "http://www.xml-cml.org/convention/");
  cml.put("<xmlattr>.convention", "convention:molecular");
}

CMLWriter::~CMLWriter() { write(); }

void CMLWriter::write() const {
  PRECONDITION(p_ostream, "no output stream");
  boost::property_tree::write_xml(
      *p_ostream, tree,
      boost::property_tree::xml_writer_make_settings<std::string>(' ', 2));
}

void CMLWriter::add_molecule(const ROMol &mol, int confId) {
  auto &molecule_node = tree.add("cml.molecule", "");
  molecule_node.put(
      "<xmlattr>.id",
      boost::format{"%1%%2%"} % molecule_id_prefix % num_written_mols++);

  RWMol rwmol{mol};
  put_atomArray(molecule_node, rwmol, confId);
  put_bondArray(molecule_node);
}

void CMLWriter::put_atomArray(boost::property_tree::ptree &molecule_node,
                              const RWMol &rwmol, int confId) {
  const Conformer *conf =
      rwmol.getNumConformers() ? &rwmol.getConformer(confId) : nullptr;

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
}

void CMLWriter::put_bondArray(boost::property_tree::ptree &molecule_node) {
  unsigned bond_id = 0u;
  auto &bondArray = molecule_node.put("bondArray", "");
#if 0
  for (auto atom_itr = rwmol.beginAtoms(), atom_itr_end = rwmol.endAtoms();
       atom_itr != atom_itr_end; ++atom_itr) {
    const auto &atom = *atom_itr;
    PRECONDITION(atom, "bad atom");
    const auto src = atom->getIdx();
    for (auto bond_itrs = rwmol.getAtomBonds(atom);
         bond_itrs.first != bond_itrs.second; ++bond_itrs.first) {
      auto *bptr = rwmol[*bond_itrs.first];
      auto *nptr = bptr->getOtherAtom(atom);
      const auto dst = nptr->getIdx();
      if (dst < src) {
        continue;
      }

      auto &bond = bondArray.add("bond", "");
      bond.put("<xmlattr>.atomRefs2",
               boost::format{"%1%%2% %1%%3%"} % atom_id_prefix % src % dst);

      bond.put("<xmlattr>.id",
               boost::format{"%1%%2%"} % bond_id_prefix % bond_id++);
    }
  }
#endif
}
}  // namespace RDKit
