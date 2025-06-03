#include "CJson.h"
#include <iostream>
#include "GraphMol/FileParsers/FileParsers.h"

int main() {
  const RDKit::v2::FileParsers::MolFileParserParams params{.sanitize = false,
                                                           .removeHs = false};

  auto mol = RDKit::v2::FileParsers::MolFromMolBlock(
#if 1
      R"(
 OpenBabel06022510293D

 10  9  0  0  0  0  0  0  0  0999 V2000
    0.9704   -0.0287    0.0789 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4707    1.3578   -0.3862 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4536    2.1227   -0.6413 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7670    1.5178   -0.4413 O   0  5  0  0  0  0  0  0  0  0  0  0
    2.4790   -0.0695    0.0926 N   0  3  0  0  0  0  0  0  0  0  0  0
    0.6419   -0.8145   -0.6058 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6419   -0.2425    1.0991 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.7711    0.2306   -0.8491 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.9076   -0.9434    0.3864 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.7708    0.7385    0.6622 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  5  1  0  0  0  0
  1  6  1  0  0  0  0
  1  7  1  0  0  0  0
  2  3  2  0  0  0  0
  2  4  1  0  0  0  0
  5  8  1  0  0  0  0
  5  9  1  0  0  0  0
  5 10  1  0  0  0  0
M  CHG  2   4  -1   5   1
M  END
)"
#endif
      ,
      params);

  // RDKit::MolToCJSONBlock(std::cout, *mol);
  // std::cout << std::endl;
  RDKit::MolToCJSONFile(*mol, "zzz.json");
}
