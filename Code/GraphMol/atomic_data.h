//
//  Copyright (C) 2001-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*! \file atomic_data.h

  \brief No user-serviceable parts inside

  This stuff is used by the PeriodicTable interface

*/
#include <RDGeneral/export.h>
#ifndef __RD_ATOMIC_DATA_H
#define __RD_ATOMIC_DATA_H

#include <RDGeneral/types.h>
#include <map>

namespace RDKit {
RDKIT_GRAPHMOL_EXPORT extern const std::string periodicTableAtomData;
RDKIT_GRAPHMOL_EXPORT extern const std::string isotopesAtomData[];
RDKIT_GRAPHMOL_EXPORT extern const std::vector<std::string> elementNames;
namespace constants {
RDKIT_GRAPHMOL_EXPORT extern const double electronMass;
}
class RDKIT_GRAPHMOL_EXPORT atomicData {
 public:
  atomicData(const std::string &dataLine);
  ~atomicData() = default;

  [[nodiscard]] int AtomicNum() const { return anum; }

  [[nodiscard]] int DefaultValence() const { return valence.front(); }

  [[nodiscard]] int NumValence() const {
    return static_cast<int>(valence.size());
  }

  [[nodiscard]] const INT_VECT &ValenceList() const { return valence; }

  [[nodiscard]] double Mass() const { return mass; }

  [[nodiscard]] std::string Symbol() const { return symb; }

  [[nodiscard]] std::string Name() const { return name; }

  [[nodiscard]] double Rcov() const { return rCov; }

  [[nodiscard]] double Rb0() const { return rB0; }

  [[nodiscard]] double Rvdw() const { return rVdw; }

  [[nodiscard]] int NumOuterShellElec() const { return nVal; }

  [[nodiscard]] int MostCommonIsotope() const { return commonIsotope; }

  [[nodiscard]] double MostCommonIsotopeMass() const {
    return commonIsotopeMass;
  }

  [[nodiscard]] unsigned int Row() const { return row; }
  // maps isotope number -> mass
  std::map<unsigned int, std::pair<double, double>>
      d_isotopeInfoMap;  // available isotopes
 private:
  int anum;                // atomic number
  std::string symb;        // atomic symbol
  std::string name;        // atomic name
  double rCov, rB0, rVdw;  // radii
  INT_VECT valence;        // list of all valences, the first one is the default
  // valence, -1 at the end signifies that any upper valence
  // is tolerated
  double mass;               // atomic mass
  int nVal;                  // number of outer shell electrons
  int commonIsotope;         // most common isotope
  double commonIsotopeMass;  // most common isotopic mass
  unsigned int row;          // row in the periodic table
};
};  // namespace RDKit
#endif
