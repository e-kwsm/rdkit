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
#include "CMLReader.h"

#include <boost/format.hpp>

#ifndef __clang__version__    // FIXME
#pragma GCC diagnostic pop    // FIXME
#else                         // FIXME
#pragma clang diagnostic pop  // FIXME
#endif                        // FIXME

using namespace std::literals::string_literals;
using namespace RDKit;
using namespace RDKit::v2::FileParsers;

namespace {}  // namespace

SCENARIO("CML Reader", "[CML][reader]") {
  using namespace RDKit::v2::FileParsers;
  WHEN("unreadable file is passed") { REQUIRE_THROWS(CMLSupplier{"/a.cml"}); }

  WHEN("XML is malformed") {
    {
      std::stringbuf buf{R"(<?xml version="1.0"?><cml/><cml/>)"s};
      std::unique_ptr<std::istream> pis = std::make_unique<std::istream>(&buf);
      REQUIRE_THROWS_AS(CMLSupplier{std::move(pis)}, RDKit::FileParseException);
    }
    {
      std::stringbuf buf{R"(<?xml version="1.0"?><cml>)"s};
      std::unique_ptr<std::istream> pis = std::make_unique<std::istream>(&buf);
      REQUIRE_THROWS_AS(CMLSupplier{std::move(pis)}, RDKit::FileParseException);
    }
  }

  WHEN("//molecule/atomArray has no atom elements") {
    std::stringbuf buf{R"(<?xml version="1.0"?>
<cml>
  <molecule>
    <atomArray/>
  </molecule>
</cml>
)"s};
    std::unique_ptr<std::istream> pis = std::make_unique<std::istream>(&buf);
    CMLSupplier supp{std::move(pis)};
    REQUIRE_THROWS_AS(supp.next(), RDKit::FileParseException);

    // auto cml = cml_root(MoleculeNode{}.str(cml_atomArray("")));
    // REQUIRE_THROWS_WITH(
    //     MolFromCMLBlock(cml),
    //     Matches(".*/molecule[^/]*/atomArray has no atom elements"));
  }
}
