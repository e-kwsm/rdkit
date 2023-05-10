#include <stdexcept>

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

#include "RDGeneral/test.h"

#define CATCH_CONFIG_MAIN  // FIXME

#include "catch.hpp"

#include <GraphMol/FileParsers/FileParsers.h>

#include <boost/format.hpp>

#ifndef __clang__version__    // FIXME
#pragma GCC diagnostic pop    // FIXME
#else                         // FIXME
#pragma clang diagnostic pop  // FIXME
#endif                        // FIXME

using namespace RDKit;

#define LOCATION                                           \
  do {                                                     \
    std::cerr << __FILE__ << ':' << __LINE__ << std::endl; \
  } while (0)  // FIXME

// FIXME

SCENARIO("CML Reader", "[CML][reader]") {
  // RDLog::InitLogs();  // FIXME
  using namespace std::literals::string_literals;
  using Catch::Matches;

  // Use functor since reference to lambda, a local variable, is not allowed
  struct put_unless_none {
    void operator()(std::ostream& os, std::string key,
                    const boost::optional<std::string>& val) {
      if (val) {
        os << boost::format{" %1%='%2%'"} % key % *val;
      }
    }
  };

  auto cml_root = [](const std::string& body) {
    return "<cml"
           " xmlns='http://www.xml-cml.org/schema'"
           " xmlns:convention='http://www.xml-cml.org/convention/'"
           " convention='convention:molecular'>" +
           body + "</cml>";
  };

  struct molecule_node {
#if 0 && __cplusplus < 201402L
//    // XXX C++11 requires the following constructors
//    molecule_node() = default;
//    molecule_node(std::string id, std::string formalCharge,
//                  std::string spinMultiplicity)
//        : id{std::move(id)},
//          formalCharge{std::move(formalCharge)},
//          spinMultiplicity{std::move(spinMultiplicity)} {}
#endif

    boost::optional<std::string> id{"m0"};
    boost::optional<std::string> formalCharge{"0"};
    boost::optional<std::string> spinMultiplicity{"1"};

    std::string str(std::string body = "") const {
      boost::format fmt{" %1%='%2%'"};
      std::ostringstream ss;
      ss << "<molecule";
      put_unless_none{}(ss, "id", id);
      put_unless_none{}(ss, "formalCharge", formalCharge);
      put_unless_none{}(ss, "spinMultiplicity", spinMultiplicity);
      ss << ">" << body << "</molecule>";
      return ss.str();
    }
  };

  auto cml_atomArray = [](std::string body) {
    return "<atomArray>" + body + "</atomArray>";
  };

  struct atom_node {
#if 0 && __cplusplus < 201402L
    // XXX C++11 requires the following constructors
    atom_node() = default;
    atom_node(std::string id, std::string elementType, std::string formalCharge,
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

    boost::optional<std::string> id{"a0"};
    boost::optional<std::string> elementType{"He"};
    boost::optional<std::string> formalCharge{"0"};
    boost::optional<std::string> spinMultiplicity{"1"};
    boost::optional<std::string> hydrogenCount{"0"};
    boost::optional<std::string> x3{"0.0"}, y3{"0.0"}, z3{"0.0"};

    auto str() const {
      std::ostringstream ss;
      ss << "<atom";
      put_unless_none{}(ss, "elementType", elementType);
      put_unless_none{}(ss, "id", id);
      put_unless_none{}(ss, "formalCharge", formalCharge);
      put_unless_none{}(ss, "spinMultiplicity", spinMultiplicity);
      put_unless_none{}(ss, "hydrogenCount", hydrogenCount);
      put_unless_none{}(ss, "x3", x3);
      put_unless_none{}(ss, "y3", y3);
      put_unless_none{}(ss, "z3", z3);
      ss << "/>";
      return ss.str();
    }
  };

  auto cml_bondArray = [](std::string body) {
    return "<bondArray>" + body + "</bondArray>";
  };

  struct bond_node {
#if 0 && __cplusplus < 201402L
    // XXX C++11 requires the following constructors
    bond_node(std::string id, std::string atomRefs2, std::string order)
        : id{std::move(id)},
          atomRefs2{std::move(atomRefs2)},
          order{std::move(order)} {}
#endif

    boost::optional<std::string> id{"b0"};
    boost::optional<std::string> atomRefs2{"a0 a1"};
    boost::optional<std::string> order{"1"};
    // TODO
    // boost::optional<std::string> bondStereo{boost::none};

    auto str() const {
      std::ostringstream ss;
      ss << "<bond";
      put_unless_none{}(ss, "id", id);
      put_unless_none{}(ss, "atomRefs2", atomRefs2);
      put_unless_none{}(ss, "order", order);
      ss << "/>";
      return ss.str();
    }
  };

  //

  WHEN("multiple root nodes exist") {
    REQUIRE_THROWS(CMLBlockToMol(cml_root("") + cml_root("")));
  }

  WHEN("//molecule/@formalCharge") {
    AND_WHEN("is not numeric") {
      molecule_node m;
      m.formalCharge = "abc";
      auto cml = cml_root(m.str());
      auto msg =
          boost::format{
              ".*/molecule[^/]*/@formalCharge "
              R"CML(\(= "%s"\) is not convertible to numeric)CML"} %
          *m.formalCharge;
      REQUIRE_THROWS_WITH(CMLBlockToMol(cml), Matches(msg.str()));
    }
  }

  WHEN("//molecule/@spinMultiplicity") {
    AND_WHEN("is not numeric") {
      molecule_node m;
      m.spinMultiplicity = "abc";
      auto cml = cml_root(m.str());
      auto msg =
          boost::format{
              ".*/molecule[^/]*/@spinMultiplicity "
              R"CML(\(= "%s"\) is not convertible to numeric)CML"} %
          *m.spinMultiplicity;
      REQUIRE_THROWS_WITH(CMLBlockToMol(cml), Matches(msg.str()));
    }
    AND_WHEN("is zero") {
      molecule_node m;
      m.spinMultiplicity = "0";
      auto cml = cml_root(m.str());
      REQUIRE_THROWS_WITH(
          CMLBlockToMol(cml),
          Matches(".*/molecule[^/]*/@spinMultiplicity is zero"));
    }
  }

  WHEN("//molecule/atomArray") {
    AND_WHEN("has no atom elements") {
      auto cml = cml_root(molecule_node{}.str(cml_atomArray("")));
      REQUIRE_THROWS_WITH(
          CMLBlockToMol(cml),
          Matches(".*/molecule[^/]*/atomArray has no atom elements"));
    }
  }

  WHEN("//molecule/atomArray/atom/@id") {
    AND_WHEN("is missing") {
      atom_node a;
      a.id = boost::none;
      auto cml = cml_root(molecule_node{}.str(cml_atomArray(a.str())));
      REQUIRE_THROWS_WITH(CMLBlockToMol(cml),
                          Matches(".+/atom/@id is missing"));
    }

    AND_WHEN("is invalid") {
      atom_node a;
      a.id = "a0 a1";
      auto cml = cml_root(molecule_node{}.str(cml_atomArray(a.str())));
      auto msg =
          boost::format{R"CML(.+/atom/@id .*\(= "%s"\) is invalid)CML"} % *a.id;
      REQUIRE_THROWS_WITH(CMLBlockToMol(cml), Matches(msg.str()));
    }

    AND_WHEN("is not unique") {
      atom_node a;
      auto cml =
          cml_root(molecule_node{}.str(cml_atomArray(a.str() + a.str())));
      auto msg =
          boost::format{R"CML(.+/atom/@id .*\(= "%s"\) is not unique)CML"} %
          *a.id;
      REQUIRE_THROWS_WITH(CMLBlockToMol(cml), Matches(msg.str()));
    }
  }

  WHEN("//molecule/atomArray/atom/@elementType") {
    AND_WHEN("is missing") {
      atom_node a;
      a.elementType = boost::none;
      auto cml = cml_root(molecule_node{}.str(cml_atomArray(a.str())));
      REQUIRE_THROWS_WITH(CMLBlockToMol(cml),
                          Matches(".+/atom[^/]*/@elementType is missing"));
    }
  }

  WHEN("//molecule/atomArray/atom/@formalCharge") {
    AND_WHEN("disagrees with //molecule/@formalCharge") {
      molecule_node m;
      m.formalCharge = "-1";
      atom_node a;
      a.formalCharge = "0";
      auto cml = cml_root(m.str(cml_atomArray(a.str())));
      auto msg =
          boost::format{
              ".*/molecule[^/]*/@formalCharge "
              R"CML(\(= %d\))CML"
              " is not equal to sum of .+/atom/@formalCharge "
              R"CML(\(= %d\))CML"} %
          *m.formalCharge % *a.formalCharge;
      REQUIRE_THROWS_WITH(CMLBlockToMol(cml), Matches(msg.str()));
    }
  }

  WHEN("//molecule/atomArray/atom/@spinMultiplicity") {
    AND_WHEN("is zero") {
      atom_node a;
      a.spinMultiplicity = "0";
      auto cml = cml_root(molecule_node{}.str(cml_atomArray(a.str())));
      REQUIRE_THROWS_WITH(CMLBlockToMol(cml),
                          Matches(".+/atom[^/]*/@spinMultiplicity is zero"));
    }
#if 0
    AND_WHEN("disagrees with //molecule/@spinMultiplicity") {
      molecule_node m;
      m.spinMultiplicity = "3";
      atom_node a;
      a.spinMultiplicity = "1";
      auto cml = cml_root(m.str(cml_atomArray(a.str())));
      REQUIRE_THROWS_WITH(
          CMLBlockToMol(cml),
          Matches(".*/molecule[^/]*/@spinMultiplicity .+ is not equal to .* "
                  ".+/atom/@spinMultiplicity.+"));
    }
#endif
  }

  WHEN("//molecule/atomArray/atom/@x3") {
    AND_WHEN("is missing") {
      atom_node a;
      a.x3 = boost::none;
      auto cml = cml_root(molecule_node{}.str(cml_atomArray(a.str())));
      REQUIRE_THROWS_WITH(
          CMLBlockToMol(cml),
          Matches(
              ".+/atom[^/]* does not have all of x3, y3 and z3 attributes"));
    }
  }

  WHEN("//molecule/bondArray") {
    AND_WHEN("has no bond elements") {
      auto cml = cml_root(molecule_node{}.str(cml_atomArray(atom_node{}.str()) +
                                              cml_bondArray("")));
      REQUIRE_THROWS_WITH(
          CMLBlockToMol(cml),
          Matches(".+/molecule[^/]*/bondArray has no bond elements"));
    }
  }

  // hydrogen cyanide
  const atom_node a0{"a0"s, "H"s, "0"s, "1"s, "0"s, "-1.06"s, "0.0"s, "0.0"s};
  const atom_node a1{"a1"s, "C"s, "0"s, "1"s, "1"s, "0.00"s, "0.0"s, "0.0"s};
  const atom_node a2{"a2"s, "N"s, "0"s, "1"s, "0"s, "1.16"s, "0.0"s, "0.0"s};
  const auto atomArray = cml_atomArray(a0.str() + a1.str() + a2.str());
  const bond_node b0{"b0"s, "a0 a1"s, "1"s};
  const bond_node b1{"b1"s, "a1 a2"s, "3"s};

  WHEN("//molecule/bondArray/bond/@id") {
    AND_WHEN("is not unique") {
      auto b = b0;
      b.id = b1.id;
      const auto cml = cml_root(
          molecule_node{}.str(atomArray + cml_bondArray(b.str() + b1.str())));
      REQUIRE_THROWS_WITH(CMLBlockToMol(cml), Matches(".+/bond/@id .*" + *b.id +
                                                      ".* is not unique"));
    }
  }

  WHEN("//molecule/bondArray/bond/@atomRefs2") {
    AND_WHEN("is missing") {
      auto b = b0;
      b.atomRefs2 = boost::none;
      const auto cml = cml_root(
          molecule_node{}.str(atomArray + cml_bondArray(b.str() + b1.str())));
      REQUIRE_THROWS_WITH(CMLBlockToMol(cml),
                          Matches(".+/bond[^/]*/@atomRefs2 is missing"));
    }

    AND_WHEN("has a single atom id") {
      auto b = b0;
      b.atomRefs2 = "a1";
      const auto cml = cml_root(
          molecule_node{}.str(atomArray + cml_bondArray(b.str() + b1.str())));
      REQUIRE_THROWS_WITH(
          CMLBlockToMol(cml),
          Matches(".+/bond[^/]*/@atomRefs2 .*" + *b.atomRefs2 +
                  ".* does not have two ids separated by space"));
    }

    AND_WHEN("has too many atom ids") {
      auto b = b0;
      b.atomRefs2 = "a0 a1 a2";
      const auto cml = cml_root(
          molecule_node{}.str(atomArray + cml_bondArray(b.str() + b1.str())));
      REQUIRE_THROWS_WITH(CMLBlockToMol(cml),
                          Matches(".+/bond[^/]*/@atomRefs2 .*" + *b.atomRefs2 +
                                  ".* has three or more ids"));
    }

    AND_WHEN("is self-bond") {
      auto b = b0;
      b.atomRefs2 = "a0 a0";
      const auto cml = cml_root(
          molecule_node{}.str(atomArray + cml_bondArray(b.str() + b1.str())));
      REQUIRE_THROWS_WITH(CMLBlockToMol(cml),
                          Matches(".+/bond[^/]*/@atomRefs2 .*" + *b.atomRefs2 +
                                  ".* is self-bond"));
    }

    AND_WHEN("refers to non-existing //atom/@id") {
      auto a = "a-1";
      auto b = b0;
      b.atomRefs2 = *a0.id + " " + a;
      const auto cml = cml_root(
          molecule_node{}.str(atomArray + cml_bondArray(b.str() + b1.str())));
      REQUIRE_THROWS_WITH(CMLBlockToMol(cml),
                          Matches(".+/bond[^/]*/@atomRefs2 .*" + *b.atomRefs2 +
                                  ".* refers to non-existing .+" + a + ".*"));
    }
  }

  WHEN("//molecule/bondArray/bond/@order") {
    AND_WHEN("is missing") {
      auto b = b0;
      b.order = boost::none;
      const auto cml = cml_root(
          molecule_node{}.str(atomArray + cml_bondArray(b.str() + b1.str())));
      REQUIRE_THROWS_WITH(CMLBlockToMol(cml),
                          Matches(".+/bond[^/]*/@order is missing"));
    }

#if 0
    AND_WHEN("is invalid") {
      bond_attributes attr;
      attr.order = "-1";
      const auto cml = cml_root(cml_molecule(
          cml_atomArray(cml_atom({})) + cml_bondArray(cml_bond(attr)), {}));
      REQUIRE_NOTHROW(CMLBlockToMol(cml));
    }
#endif
  }
}
