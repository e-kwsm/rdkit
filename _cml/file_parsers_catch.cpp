#ifndef __clang__version__                     // FIXME
#pragma GCC diagnostic push                    // FIXME
#pragma GCC diagnostic ignored "-Wall"         // FIXME
#pragma GCC diagnostic ignored "-Wextra"       // FIXME
#pragma GCC diagnostic ignored "-Wpedantic"    // FIXME
#else                                          // FIXME
#pragma clang diagnostic push                  // FIXME
#pragma clang diagnostic ignored "-Wall"       // FIXME
#pragma clang diagnostic ignored "-Wextra"     // FIXME
#pragma clang diagnostic ignored "-Wpedantic"  // FIXME
#endif                                         // FIXME

#include <optional>
#include <stdexcept>

#include "RDGeneral/test.h"
#include <catch2/catch_all.hpp>
#include <GraphMol/FileParsers/FileParsers.h>

#include <sstream>
#include <boost/format.hpp>

#ifndef __clang__version__    // FIXME
#pragma GCC diagnostic pop    // FIXME
#else                         // FIXME
#pragma clang diagnostic pop  // FIXME
#endif                        // FIXME

using namespace RDKit;
using RDKit::v2::FileParsers::MolFromCMLBlock;

#define LOCATION                                           \
  do {                                                     \
    std::cerr << __FILE__ << ':' << __LINE__ << std::endl; \
  } while (0)  // FIXME

// FIXME

namespace {
std::string put_attribute_unless_none(const std::string& key,
                                      const std::optional<std::string>& val) {
  return val ? (boost::format{" %1%='%2%'"} % key % *val).str() : "";
}

auto cml_root(const std::string& body) {
  return "<cml"
         " xmlns='http://www.xml-cml.org/schema'"
         " xmlns:convention='http://www.xml-cml.org/convention/'"
         " convention='convention:molecular'>" +
         body + "</cml>";
}

struct MoleculeNode {
#if 0 && __cplusplus < 201402L
  // XXX C++11 requires the following constructors
  MoleculeNode() = default;
  MoleculeNode(std::string id, std::string formalCharge,
               std::string spinMultiplicity)
      : id{std::move(id)},
        formalCharge{std::move(formalCharge)},
        spinMultiplicity{std::move(spinMultiplicity)} {}
#endif

  std::optional<std::string> id{"m0"};
  std::optional<std::string> formalCharge{"0"};
  std::optional<std::string> spinMultiplicity{"1"};

  std::string str(const std::string& body = "") const {
    boost::format fmt{" %1%='%2%'"};
    std::ostringstream ss;
    ss << "<molecule";
    ss << put_attribute_unless_none("id", id);
    ss << put_attribute_unless_none("formalCharge", formalCharge);
    ss << put_attribute_unless_none("spinMultiplicity", spinMultiplicity);
    ss << ">" << body << "</molecule>";
    return ss.str();
  }
};

auto cml_atomArray(const std::string& body) {
  return "<atomArray>" + body + "</atomArray>";
}

struct AtomNode {
#if 0 && __cplusplus < 201402L
  // XXX C++11 requires the following constructors
  AtomNode() = default;
  AtomNode(std::string id, std::string elementType, std::string formalCharge,
           std::string spinMultiplicity, std::string hydrogenCount,
           std::string x3, std::string y3, std::string z3)
      : id{std::move(id)},
        elementType{std::move(elementType)},
        formalCharge{std::move(formalCharge)},
        spinMultiplicity{std::move(spinMultiplicity)},
        hydrogenCount{std::move(hydrogenCount)},
        x3{std::move(x3)},
        y3{std::move(y3)},
        z3{std::move(z3)} {}
#endif

  std::optional<std::string> id{"a0"};
  std::optional<std::string> elementType{"He"};
  std::optional<std::string> formalCharge{"0"};
  std::optional<std::string> spinMultiplicity{"1"};
  std::optional<std::string> hydrogenCount{"0"};
  std::optional<std::string> x3{"0.0"}, y3{"0.0"}, z3{"0.0"};

  std::string str() const {
    std::ostringstream ss;
    ss << "<atom";
    ss << put_attribute_unless_none("elementType", elementType);
    ss << put_attribute_unless_none("id", id);
    ss << put_attribute_unless_none("formalCharge", formalCharge);
    ss << put_attribute_unless_none("spinMultiplicity", spinMultiplicity);
    ss << put_attribute_unless_none("hydrogenCount", hydrogenCount);
    ss << put_attribute_unless_none("x3", x3);
    ss << put_attribute_unless_none("y3", y3);
    ss << put_attribute_unless_none("z3", z3);
    ss << "/>";
    return ss.str();
  }
};

auto cml_bondArray(const std::string& body) {
  return "<bondArray>" + body + "</bondArray>";
}

struct BondNode {
#if 0 && __cplusplus < 201402L
  // XXX C++11 requires the following constructors
  BondNode(std::string id, std::string atomRefs2, std::string order)
      : id{std::move(id)},
        atomRefs2{std::move(atomRefs2)},
        order{std::move(order)} {}
#endif

  std::optional<std::string> id{"b0"};
  std::optional<std::string> atomRefs2{"a0 a1"};
  std::optional<std::string> order{"1"};
  // TODO
  // std::optional<std::string> bondStereo{std::nullopt};

  std::string str() const {
    std::ostringstream ss;
    ss << "<bond";
    ss << put_attribute_unless_none("id", id);
    ss << put_attribute_unless_none("atomRefs2", atomRefs2);
    ss << put_attribute_unless_none("order", order);
    ss << "/>";
    return ss.str();
  }
};
}  // namespace

SCENARIO("CML Reader", "[CML][reader]") {
  // RDLog::InitLogs();  // FIXME
  using namespace std::literals::string_literals;
  using Catch::Matchers::Matches;

  WHEN("multiple root nodes exist") {
    REQUIRE_THROWS(MolFromCMLBlock(cml_root("") + cml_root("")));
  }

  WHEN("//molecule/@spinMultiplicity") {
    AND_WHEN("is not numeric") {
      MoleculeNode m;
      m.spinMultiplicity = "abc";
      const auto cml = cml_root(m.str());
      REQUIRE_THROWS_WITH(
          MolFromCMLBlock(cml),
          Matches(
              ".*/molecule[^/]*/@spinMultiplicity .+ is not convertible.*"));
    }
    AND_WHEN("is zero") {
      MoleculeNode m;
      m.spinMultiplicity = "0";
      const auto cml = cml_root(m.str());
      REQUIRE_THROWS_WITH(
          MolFromCMLBlock(cml),
          Matches(".*/molecule[^/]*/@spinMultiplicity is zero"));
    }
  }

  WHEN("//molecule/atomArray") {
    AND_WHEN("has no atom elements") {
      const auto cml = cml_root(MoleculeNode{}.str(cml_atomArray("")));
      REQUIRE_THROWS_WITH(MolFromCMLBlock(cml),
                          Matches(".*/molecule[^/]*/atomArray has no atom"));
    }
  }

  WHEN("//molecule/atomArray/atom/@id") {
    AND_WHEN("is missing") {
      AtomNode a;
      a.id = std::nullopt;
      const auto cml = cml_root(MoleculeNode{}.str(cml_atomArray(a.str())));
      REQUIRE_THROWS_WITH(MolFromCMLBlock(cml),
                          Matches(".+/atom/@id is missing"));
    }

    AND_WHEN("is invalid") {
      AtomNode a;
      a.id = "a0 a1";
      const auto cml = cml_root(MoleculeNode{}.str(cml_atomArray(a.str())));
      REQUIRE_THROWS_WITH(MolFromCMLBlock(cml),
                          Matches(".+/atom/@id .*" + *a.id + ".* is invalid"));
    }

    AND_WHEN("is not unique") {
      AtomNode a;
      const auto cml =
          cml_root(MoleculeNode{}.str(cml_atomArray(a.str() + a.str())));
      REQUIRE_THROWS_WITH(
          MolFromCMLBlock(cml),
          Matches(".+/atom/@id .*" + *a.id + ".* is not unique"));
    }
  }

  WHEN("//molecule/atomArray/atom/@elementType") {
    AND_WHEN("is missing") {
      AtomNode a;
      a.elementType = std::nullopt;
      const auto cml = cml_root(MoleculeNode{}.str(cml_atomArray(a.str())));
      REQUIRE_THROWS_WITH(MolFromCMLBlock(cml),
                          Matches(".+/atom[^/]*/@elementType is missing"));
    }
  }

  WHEN("//molecule/atomArray/atom/@formalCharge") {
    AND_WHEN("disagrees with //molecule/@formalCharge") {
      MoleculeNode m;
      m.formalCharge = "-1";
      AtomNode a;
      a.formalCharge = "0";
      const auto cml = cml_root(m.str(cml_atomArray(a.str())));
      REQUIRE_THROWS_WITH(
          MolFromCMLBlock(cml),
          Matches(".*/molecule[^/]*/@formalCharge .+ is not equal to .* "
                  ".+/atom/@formalCharge.+"));
    }
  }

  WHEN("//molecule/atomArray/atom/@spinMultiplicity") {
    AND_WHEN("is zero") {
      AtomNode a;
      a.spinMultiplicity = "0";
      const auto cml = cml_root(MoleculeNode{}.str(cml_atomArray(a.str())));
      REQUIRE_THROWS_WITH(MolFromCMLBlock(cml),
                          Matches(".+/atom[^/]*/@spinMultiplicity is zero"));
    }
    AND_WHEN("disagrees with //molecule/@spinMultiplicity") {
      MoleculeNode m;
      m.spinMultiplicity = "3";
      AtomNode a;
      a.spinMultiplicity = "1";
      const auto cml = cml_root(m.str(cml_atomArray(a.str())));
      REQUIRE_THROWS_WITH(
          MolFromCMLBlock(cml),
          Matches(".*/molecule[^/]*/@spinMultiplicity .+ is not equal to .* "
                  ".+/atom/@spinMultiplicity.+"));
    }
  }

  WHEN("//molecule/atomArray/atom/@x3") {
    AND_WHEN("is missing") {
      AtomNode a;
      a.x3 = std::nullopt;
      const auto cml = cml_root(MoleculeNode{}.str(cml_atomArray(a.str())));
      REQUIRE_THROWS_WITH(
          MolFromCMLBlock(cml),
          Matches(
              ".+/atom[^/]* does not have all of x3, y3 and z3 attributes"));
    }
  }

  WHEN("//molecule/bondArray") {
    AND_WHEN("has no bond elements") {
      const auto cml = cml_root(MoleculeNode{}.str(
          cml_atomArray(AtomNode{}.str()) + cml_bondArray("")));
      REQUIRE_THROWS_WITH(MolFromCMLBlock(cml),
                          Matches(".+/molecule[^/]*/bondArray has no bond"));
    }
  }

  // hydrogen cyanide
  const AtomNode a0{"a0"s, "H"s, "0"s, "1"s, "0"s, "-1.06"s, "0.0"s, "0.0"s};
  const AtomNode a1{"a1"s, "C"s, "0"s, "1"s, "1"s, "0.00"s, "0.0"s, "0.0"s};
  const AtomNode a2{"a2"s, "N"s, "0"s, "1"s, "0"s, "1.16"s, "0.0"s, "0.0"s};
  const auto atomArray = cml_atomArray(a0.str() + a1.str() + a2.str());
  const BondNode b0{"b0"s, "a0 a1"s, "1"s};
  const BondNode b1{"b1"s, "a1 a2"s, "3"s};

  WHEN("//molecule/bondArray/bond/@id") {
    AND_WHEN("is not unique") {
      auto b = b0;
      b.id = b1.id;
      const auto cml = cml_root(
          MoleculeNode{}.str(atomArray + cml_bondArray(b.str() + b1.str())));
      REQUIRE_THROWS_WITH(
          MolFromCMLBlock(cml),
          Matches(".+/bond/@id .*" + *b.id + ".* is not unique"));
    }
  }

  WHEN("//molecule/bondArray/bond/@atomRefs2") {
    AND_WHEN("is missing") {
      auto b = b0;
      b.atomRefs2 = std::nullopt;
      const auto cml = cml_root(
          MoleculeNode{}.str(atomArray + cml_bondArray(b.str() + b1.str())));
      REQUIRE_THROWS_WITH(MolFromCMLBlock(cml),
                          Matches(".+/bond[^/]*/@atomRefs2 is missing"));
    }

    AND_WHEN("has a single atom id") {
      auto b = b0;
      b.atomRefs2 = "a1";
      const auto cml = cml_root(
          MoleculeNode{}.str(atomArray + cml_bondArray(b.str() + b1.str())));
      REQUIRE_THROWS_WITH(
          MolFromCMLBlock(cml),
          Matches(".+/bond[^/]*/@atomRefs2 .*" + *b.atomRefs2 +
                  ".* does not have two ids separated by space"));
    }

    AND_WHEN("has too many atom ids") {
      auto b = b0;
      b.atomRefs2 = "a0 a1 a2";
      const auto cml = cml_root(
          MoleculeNode{}.str(atomArray + cml_bondArray(b.str() + b1.str())));
      REQUIRE_THROWS_WITH(MolFromCMLBlock(cml),
                          Matches(".+/bond[^/]*/@atomRefs2 .*" + *b.atomRefs2 +
                                  ".* has three or more ids"));
    }

    AND_WHEN("is self-bond") {
      auto b = b0;
      b.atomRefs2 = "a0 a0";
      const auto cml = cml_root(
          MoleculeNode{}.str(atomArray + cml_bondArray(b.str() + b1.str())));
      REQUIRE_THROWS_WITH(MolFromCMLBlock(cml),
                          Matches(".+/bond[^/]*/@atomRefs2 .*" + *b.atomRefs2 +
                                  ".* is self-bond"));
    }

    AND_WHEN("refers to non-existing //atom/@id") {
      auto a = "a-1";
      auto b = b0;
      b.atomRefs2 = *a0.id + " " + a;
      const auto cml = cml_root(
          MoleculeNode{}.str(atomArray + cml_bondArray(b.str() + b1.str())));
      REQUIRE_THROWS_WITH(MolFromCMLBlock(cml),
                          Matches(".+/bond[^/]*/@atomRefs2 .*" + *b.atomRefs2 +
                                  ".* refers to non-existing .+" + a + ".*"));
    }
  }

  WHEN("//molecule/bondArray/bond/@order") {
    AND_WHEN("is missing") {
      auto b = b0;
      b.order = std::nullopt;
      const auto cml = cml_root(
          MoleculeNode{}.str(atomArray + cml_bondArray(b.str() + b1.str())));
      REQUIRE_THROWS_WITH(MolFromCMLBlock(cml),
                          Matches(".+/bond[^/]*/@order is missing"));
    }

#if 0
    AND_WHEN("is invalid") {
      bond_attributes attr;
      attr.order = "-1";
      const auto cml = cml_root(cml_molecule(
          cml_atomArray(cml_atom({})) + cml_bondArray(cml_bond(attr)), {}));
      REQUIRE_NOTHROW(MolFromCMLBlock(cml));
    }
#endif
  }
}
