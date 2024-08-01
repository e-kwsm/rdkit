#pragma once

#include <iosfwd>
#include <memory>
#include <string>

#include <boost/property_tree/ptree.hpp>

#include <RDGeneral/export.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>

namespace RDKit {
class RDKIT_FILEPARSERS_EXPORT CMLWriter {
 public:
  CMLWriter() = delete;
  CMLWriter(const std::string &fileName);
  CMLWriter(std::unique_ptr<std::ostream> &&p);
  CMLWriter(const CMLWriter &) = delete;
  ~CMLWriter();

  CMLWriter &operator=(const CMLWriter &) = delete;

  /// add the molecule to the CML document
  void add_molecule(const ROMol &mol, int confId = -1);
  /// write CML
  void write() const;

 private:
  void put_atomArray(boost::property_tree::ptree &molecule_node,
                     const RWMol &rwmol, const Conformer *const conformer);
  void put_bondArray(boost::property_tree::ptree &molecule_node,
                     const RWMol &rwmol);
  static std::string bond_order(const Bond &bond);

  std::unique_ptr<std::ostream> p_ostream;
  boost::property_tree::ptree tree;
  unsigned int num_written_mols = 0u;
  static constexpr auto molecule_id_prefix = "m";
  static constexpr auto atom_id_prefix = "a";
  static constexpr auto bond_id_prefix = "b";
};
}  // namespace RDKit
