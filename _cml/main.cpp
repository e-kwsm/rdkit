#include <istream>
#pragma GCC diagnostic push                  // FIXME
#pragma GCC diagnostic ignored "-Wall"       // FIXME
#pragma GCC diagnostic ignored "-Wextra"     // FIXME
#pragma GCC diagnostic ignored "-Wpedantic"  // FIXME

#ifdef __clang__version__                      // FIXME
#pragma clang diagnostic push                  // FIXME
#pragma clang diagnostic ignored "-Wall"       // FIXME
#pragma clang diagnostic ignored "-Wextra"     // FIXME
#pragma clang diagnostic ignored "-Wpedantic"  // FIXME
#endif                                         // FIXME

#include <iostream>

#include <GraphMol/FileParsers/FileParsers.h>
#include "CMLReader.h"

#pragma GCC diagnostic pop    // FIXME
#ifdef __clang__version__     // FIXME
#pragma clang diagnostic pop  // FIXME
#endif                        // FIXME

using namespace RDKit;
using namespace RDKit::v2::FileParsers;

namespace {
// std::istream &&to_stream(std::string &&s) {
//   std::stringbuf buf{s};
//   std::istream is{&buf};
//   return std::move(is);
// }
// std::istream &&to_stream(std::stringstream &&s) { }
}  // namespace

int main(int argc, char **argv) {
  {
    std::stringstream ss;
    ss << R"(<?xml version="1.0"?>
<cml>
  <molecule id="m1" formalCharge="0" spinMultiplicity="1">
    <atomArray>
      <atom id="a1" elementType="F" hydrogenCount="0" formalCharge="0" spinMultiplicity="1" x3="-0.678800" y3="-1.175502" z3="0.000000"/>
      <atom id="a2" elementType="C" hydrogenCount="1" formalCharge="0" spinMultiplicity="1" x3="-0.000000" y3="0.000000" z3="0.000000"/>
      <atom id="a3" elementType="C" hydrogenCount="1" formalCharge="0" spinMultiplicity="1" x3="1.328807" y3="0.000000" z3="0.000000"/>
      <atom id="a4" elementType="F" hydrogenCount="0" formalCharge="0" spinMultiplicity="1" x3="2.007607" y3="1.175502" z3="0.000000"/>
      <atom id="a5" elementType="H" hydrogenCount="0" formalCharge="0" spinMultiplicity="1" x3="-0.678800" y3="-1.175502" z3="0.000000"/>
      <atom id="a6" elementType="H" hydrogenCount="0" formalCharge="0" spinMultiplicity="1" x3="2.007607" y3="1.175502" z3="0.000000"/>
    </atomArray>
    <bondArray>
      <bond atomRefs2="a1 a2" order="1"/>
      <bond atomRefs2="a2 a3" order="2">
        <bondStereo>T</bondStereo>
      </bond>
      <bond atomRefs2="a3 a4" order="1"/>
      <bond atomRefs2="a2 a5" order="1"/>
      <bond atomRefs2="a3 a6" order="1"/>
    </bondArray>
  </molecule>
</cml>)";

    std::stringbuf buf{ss.str()};
    std::istream is{&buf};
    CMLSupplier supp{std::unique_ptr<std::istream>(&is)};
  }

#if 0
  RDLog::InitLogs();
  auto f = [](const auto& m) {
    if (m) {
      std::cout << MolToCMLBlock(*m);
      std::cerr << MolToXYZBlock(*m);
    }
  };

  if (argc > 1) {
    for (int i = 1; i < argc; i++) {
      std::cerr << argv[i] << "\n";
      try {
        f(MolFromCMLFile(argv[i]));
      } catch (std::exception& e) {
        std::cerr << e.what() << "\n";
      }
    }
    return 0;
  }

  f(MolFromCMLBlock(R"(<?xml version="1.0"?>
<cml>
  <molecule id="m1" formalCharge="0" spinMultiplicity="1">
    <atomArray>
      <atom id="a1" elementType="F" hydrogenCount="0" formalCharge="0" spinMultiplicity="1" x3="-0.678800" y3="-1.175502" z3="0.000000"/>
      <atom id="a2" elementType="C" hydrogenCount="1" formalCharge="0" spinMultiplicity="1" x3="-0.000000" y3="0.000000" z3="0.000000"/>
      <atom id="a3" elementType="C" hydrogenCount="1" formalCharge="0" spinMultiplicity="1" x3="1.328807" y3="0.000000" z3="0.000000"/>
      <atom id="a4" elementType="F" hydrogenCount="0" formalCharge="0" spinMultiplicity="1" x3="2.007607" y3="1.175502" z3="0.000000"/>
      <atom id="a5" elementType="H" hydrogenCount="0" formalCharge="0" spinMultiplicity="1" x3="-0.678800" y3="-1.175502" z3="0.000000"/>
      <atom id="a6" elementType="H" hydrogenCount="0" formalCharge="0" spinMultiplicity="1" x3="2.007607" y3="1.175502" z3="0.000000"/>
    </atomArray>
    <bondArray>
      <bond atomRefs2="a1 a2" order="1"/>
      <bond atomRefs2="a2 a3" order="2">
        <bondStereo>T</bondStereo>
      </bond>
      <bond atomRefs2="a3 a4" order="1"/>
      <bond atomRefs2="a2 a5" order="1"/>
      <bond atomRefs2="a3 a6" order="1"/>
    </bondArray>
  </molecule>
</cml>)"));

  return 0;

  f(MolFromCMLBlock(R"(<cml>
  <molecule id="m" formalCharge="0" spinMultiplicity="1">
    <atomArray>
      <atom id="a0" elementType="F" hydrogenCount="1" isotopeNumber="18"/>
      <atom id="a1" elementType="H"/>
    </atomArray>
    <bondArray>
      <bond id="b0" atomRefs2="a0 a1" order="1"/>
    </bondArray>
  </molecule>
</cml>)"));

  return 0;

  f(MolFromCMLBlock(R"(<cml>
  <molecule>
    <atomArray>
      <atom id="a0" elementType="Du" isotopeNumber="-2"/>
    </atomArray>
  </molecule>
</cml>)"));

  f(MolFromCMLBlock(R"(<cml>
  <molecule>
    <atomArray>
      <atom id="a0" elementType="Du"/>
      <atom id="a1" elementType="Du"/>
    </atomArray>
    <bondArray>
      <bond id="b0" order="XXX" atomRefs2="a0 a0"/>
    </bondArray>
  </molecule>
</cml>)"));

  f(MolFromCMLBlock(R"(<cml>
  <molecule>
    <atomArray/>
    <atomArray/>
  </molecule>
</cml>)"));

  f(MolFromCMLBlock(R"(<cml>
  <molecule formalCharge="zzzz">
    <atomArray>
      <atom/>
    </atomArray>
  </molecule>
</cml>)"));

  return 0;

  f(MolFromCMLFile("a.cml"));
  // f(MolFromCMLFile("b.cml"));
  f(MolFromCMLFile("c.cml"));

  f(MolFromCMLBlock(R"(<cml>
  <molecule>
    <atomArray/>
  </molecule>
</cml>)"));

  f(MolFromCMLBlock(R"(<cml>
  <molecule>
    <bondArray/>
  </molecule>
</cml>)"));

  MolFromCMLFile("XXX.cml");
#endif
}
