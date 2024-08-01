#include "CMLWriter.h"
#include "GraphMol/FileParsers/FileParsers.h"
#include <boost/property_tree/ptree_fwd.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <iostream>

#include <boost/format.hpp>

int main() {
  const RDKit::v2::FileParsers::MolFileParserParams params{.removeHs = false};
  RDKit::CMLWriter w{"a.cml"};

  w.add_molecule(*RDKit::v2::FileParsers::MolFromMolBlock(R"(Cl-
 OpenBabel08012407203D

  2  1  0  0  0  0  0  0  0  0999 V2000
    1.0990   -0.0544    0.0016 O   0  0  0  0  0  1  0  0  0  0  0  0
    2.0388   -0.0544    0.0016 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
)",
                                                          params));

  w.add_molecule(*RDKit::v2::FileParsers::MolFromMolBlock(R"(C2H2
  OpenBabel08012407213D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.9826   -0.0057    0.0720 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1830   -0.0057    0.0720 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0832   -0.0057    0.0720 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.2488   -0.0057    0.0720 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  3  0  0  0  0
  1  3  1  0  0  0  0
  2  4  1  0  0  0  0
M  ISO  1   2  14
M  END
)",
                                                          params));

  w.add_molecule(*RDKit::v2::FileParsers::MolFromMolBlock(R"(HCN
 OpenBabel08012407393D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.9775   -0.0570    0.0584 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.0427   -0.0570    0.0584 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2027   -0.0570    0.0584 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  3  0  0  0  0
M  END
)",
                                                          params));

  w.add_molecule(*RDKit::v2::FileParsers::MolFromMolBlock(R"(HNC
 OpenBabel08012407383D

  3  2  0  0  0  0  0  0  0  0999 V2000
    1.1033    0.0030    0.0726 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.1141    0.0030    0.0726 N   0  3  0  0  0  0  0  0  0  0  0  0
    3.2841    0.0030    0.0726 C   0  5  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  3  0  0  0  0
M  CHG  2   2   1   3  -1
M  END
)",
                                                          params));
}
