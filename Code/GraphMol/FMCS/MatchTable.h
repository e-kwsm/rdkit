//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#pragma once
#include <vector>
#include <stdexcept>

namespace RDKit {
namespace FMCS {
template <typename T>
class RDKIT_FMCS_EXPORT
    TArray2D {  // for scalar value types ! including bool with special STL
                // implementation (no reference to item - bitset used)
  size_t XSize;
  size_t YSize;
  std::vector<T> Data;

 public:
  TArray2D(size_t cy = 0, size_t cx = 0)
      : XSize(cx), YSize(cy), Data(cx * cy) {}
  size_t getXSize() const { return XSize; }
  size_t getYSize() const { return YSize; }
  bool empty() const { return Data.empty(); }
  void clear() {
    Data.clear();
    XSize = 0;
    YSize = 0;
  }
  void resize(size_t cy, size_t cx) {
    Data.resize(cx * cy);
    XSize = cx;
    YSize = cy;
  }
  void set(size_t row, size_t col, T val) { Data[row * XSize + col] = val; }
  T at(size_t row, size_t col) { return Data[row * XSize + col]; }
  T at(size_t row, size_t col) const { return Data[row * XSize + col]; }
};

typedef TArray2D<bool> MatchTable;  // row is index in QueryMolecule
}  // namespace FMCS
}  // namespace RDKit
