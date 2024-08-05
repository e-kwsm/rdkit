#include <catch2/catch_all.hpp>

#include <sstream>
#include <boost/property_tree/xml_parser.hpp>

#include "CMLReader.h"
#include "catch2/catch_test_macros.hpp"
// #include "RDGeneral/test.h"
#include <GraphMol/RWMol.h>

using namespace std::literals::string_literals;
using namespace RDKit;

namespace {
auto f(std::istream &is) {
  boost::property_tree::ptree pt;
  boost::property_tree::read_xml(is, pt);
  auto tmp = pt.get_child("molecule");
  return tmp;
}
}  // namespace

SCENARIO("CML Reader", "[CML][reader]") {
  WHEN("atomArray misses atom nodes") {
    auto block = R"(
<molecule>
  <atomArray/>
</molecule>
)"s;
    std::istringstream ss{block};
    v2::FileParsers::cml::CMLMolecule sup{f(ss)};
    REQUIRE_THROWS_AS(sup.parse(), v2::FileParsers::cml::CMLError);
  }

  WHEN("/molecule/atomArray/@id is missing") {
    auto block = R"(
<molecule>
  <atomArray>
    <atom elementType="H"/>
  </atomArray>
</molecule>
)"s;
    std::istringstream ss{block};
    v2::FileParsers::cml::CMLMolecule sup{f(ss)};
    REQUIRE_THROWS_AS(sup.parse(), v2::FileParsers::cml::CMLError);
  }

  WHEN("/molecule/atomArray/@id's are duplicated") {
    auto block = R"(
<molecule>
  <atomArray>
    <atom id="a1" elementType="H"/>
    <atom id="a1" elementType="He"/>
  </atomArray>
</molecule>
)"s;
    std::istringstream ss{block};
    v2::FileParsers::cml::CMLMolecule sup{f(ss)};
    REQUIRE_THROWS_AS(sup.parse(), v2::FileParsers::cml::CMLError);
  }

  WHEN("/molecule/atomArray/@elementType is missing") {
    auto block = R"(
<molecule>
  <atomArray>
    <atom id="a1"/>
  </atomArray>
</molecule>
)"s;
    std::istringstream ss{block};
    v2::FileParsers::cml::CMLMolecule sup{f(ss)};
    REQUIRE_THROWS_AS(sup.parse(), v2::FileParsers::cml::CMLError);
  }

#if 0
  WHEN("/molecule/atomArray/@elementType is invalid") {
    auto block = R"(
<molecule>
  <atomArray>
    <atom id="a1" elementType="ABC"/>
  </atomArray>
</molecule>
)"s;
    std::istringstream ss{block};
    v2::FileParsers::cml::CMLMolecule sup{f(ss)};
    REQUIRE_THROWS_AS(sup.parse(), Invar::Invariant);
  }
#endif

  WHEN("/molecule/atomArray/@x3 is missing") {
    for (auto a : {
             R"(<atom id="a1" elementType="H" x3="0" y3="0"/>)"s,
             R"(<atom id="a1" elementType="H" x3="0" z3="0"/>)"s,
             R"(<atom id="a1" elementType="H" y3="0" z3="0"/>)"s,
         }) {
      auto block = R"(
<molecule>
  <atomArray>
    )" + a + R"(
  </atomArray>
</molecule>
)"s;
      std::istringstream ss{block};
      v2::FileParsers::cml::CMLMolecule sup{f(ss)};
      REQUIRE_THROWS_AS(sup.parse(), v2::FileParsers::cml::CMLError);
    }
  }
}
