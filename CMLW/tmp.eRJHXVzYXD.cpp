#include "CMLWriter.h"
#include <fstream>
#include <utility>

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
  std::string name;
  rwmol.getPropIfPresent(common_properties::_Name, name);
  if (!name.empty()) {
    molecule_node.put("name", name);
  }

  put_atomArray(
      molecule_node, rwmol,
      rwmol.getNumConformers() ? &rwmol.getConformer(confId) : nullptr);
  put_bondArray(molecule_node, rwmol);
}

void CMLWriter::put_atomArray(boost::property_tree::ptree &molecule_node,
                              const ROMol &mol, const Conformer *const conf) {
  int mol_formal_charge = 0;
  unsigned mol_num_radical_electrons = 0u;

  auto &atomArray = molecule_node.put("atomArray", "");
  for (unsigned i = 0u, nAtoms = mol.getNumAtoms(); i < nAtoms; i++) {
    auto &atom_node = atomArray.add("atom", "");
    const auto &a = mol.getAtomWithIdx(i);

    atom_node.put("<xmlattr>.id", boost::format{"%1%%2%"} % atom_id_prefix % i);
    atom_node.put("<xmlattr>.elementType",
                  a->getAtomicNum() ? a->getSymbol() : "Du");

    if (conf != nullptr) {
      const auto &pos = conf->getAtomPos(i);
      boost::format xyz_fmt{"%.6f"};

      if (!conf->is3D()) {
        atom_node.put("<xmlattr>.x2", xyz_fmt % pos.x);
        atom_node.put("<xmlattr>.y2", xyz_fmt % pos.y);
      } else {
        atom_node.put("<xmlattr>.x3", xyz_fmt % pos.x);
        atom_node.put("<xmlattr>.y3", xyz_fmt % pos.y);
        atom_node.put("<xmlattr>.z3", xyz_fmt % pos.z);
      }
    }

    const auto isotope = a->getIsotope();
    if (isotope != 0u) {
      atom_node.put("<xmlattr>.isotopeNumber", isotope);
    }

    const auto charge = a->getFormalCharge();
    atom_node.put("<xmlattr>.formalCharge", charge);
    mol_formal_charge += charge;

    const auto n_rad_es = a->getNumRadicalElectrons();
    if (n_rad_es < 2u) {
      atom_node.put("<xmlattr>.spinMultiplicity", n_rad_es + 1u);
    }
    mol_num_radical_electrons += n_rad_es;
  }

  molecule_node.put("<xmlattr>.formalCharge", mol_formal_charge);

  if (mol_num_radical_electrons < 2u) {
    molecule_node.put("<xmlattr>.spinMultiplicity",
                      mol_num_radical_electrons + 1u);
  } else {
    BOOST_LOG(rdInfoLog)
        << "CMLWriter: Unable to determine molecule/@spinMultiplicity "
        << boost::format{"(%1% radical electrons)\n"} %
               mol_num_radical_electrons;
  }
}

void CMLWriter::put_bondArray(boost::property_tree::ptree &molecule_node,
                              const ROMol &mol, bool strict) {
  unsigned bond_id = 0u;
  auto &bondArray = molecule_node.put("bondArray", "");
  for (auto bond : mol.bonds()) {
    auto &bond_node = bondArray.add("bond", "");
    bond_node.put("<xmlattr>.id",
                  boost::format{"%1%%2%"} % bond_id_prefix % bond_id++);
    bond_node.put("<xmlattr>.atomRefs2",
                  boost::format{"%1%%2% %1%%3%"} % atom_id_prefix %
                      bond->getBeginAtomIdx() % bond->getEndAtomIdx());
    bond_node.put("<xmlattr>.order", bond_order(*bond, strict));

    // bond/@BondStereo if appropriate
    // http://www.xml-cml.org/convention/molecular#bondStereo-element
    auto bdir = bond->getBondDir();
    switch (bdir) {
      case Bond::BondDir::BEGINDASH:
        bond_node.put("<xmlattr>.bondStereo", "H");
        break;
      case Bond::BondDir::BEGINWEDGE:
        bond_node.put("<xmlattr>.bondStereo", "W");
        break;
      case Bond::BondDir::NONE:
        break;
      default:
        break;
    }
  }
}

std::string CMLWriter::bond_order(const Bond &bond, bool strict) {
  auto type = bond.getBondType();

  auto standardize = [=](auto val) {
    BOOST_LOG(rdWarningLog)
        << boost::format{"%1%: BondType %2% is not standard-conformant\n"} %
               __func__ % type;
    return strict ? "unknown" : val;
  };

  switch (type) {
    case Bond::SINGLE:
      return "S";
    case Bond::DOUBLE:
      return "D";
    case Bond::TRIPLE:
      return "T";
    case Bond::AROMATIC:
      return "A";

    case Bond::DATIVEONE:
      return standardize("0.5");
    case Bond::DATIVE:
      [[fallthrough]];
    case Bond::DATIVEL:
      [[fallthrough]];
    case Bond::DATIVER:
      return "S";

    case Bond::QUADRUPLE:
      return standardize("4");
    case Bond::QUINTUPLE:
      return standardize("5");
    case Bond::HEXTUPLE:
      return standardize("6");

    case Bond::ONEANDAHALF:
      return standardize("1.5");
    case Bond::TWOANDAHALF:
      return standardize("2.5");
    case Bond::THREEANDAHALF:
      return standardize("3.5");
    case Bond::FOURANDAHALF:
      return standardize("4.5");
    case Bond::FIVEANDAHALF:
      return standardize("5.5");

    case Bond::ZERO:
      return standardize("0");

    default:
      BOOST_LOG(rdWarningLog)
          << boost::format{"CMLWriter: Unsupported BondType %1%\n"} % type;
      return "unknown";
  }
}
}  // namespace RDKit
