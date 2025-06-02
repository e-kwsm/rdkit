#include "CJson.h"
#include "GraphMol/FileParsers/FileParsers.h"

int main() {
      auto mol = RDKit::v2::FileParsers::MolFromMolBlock(R"(
 OpenBabel08022400333D

  2  1  0  0  0  0  0  0  0  0999 V2000
    1.0097   -0.0332    0.0942 O   0  0  0  0  0  1  0  0  0  0  0  0
    1.9495   -0.0332    0.0942 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
)");
}
