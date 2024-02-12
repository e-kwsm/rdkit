

#if 0
[[deprecated]] bool ends_with(const std::string& str,
                              const std::string& suffix) {
  if (str.size() < suffix.size()) {
    return false;
  }
  return std::equal(suffix.rbegin(), suffix.rend(), str.rbegin());
}

namespace loglevel {
struct [[deprecated]] Error {};
}  // namespace loglevel

template <typename E>
[[deprecated]] void handler(const std::string& msg) {
  throw E{msg};
}

template <>
[[deprecated]] void handler<loglevel::Error>(const std::string& msg) {
  BOOST_LOG(rdErrorLog) << msg;
}
#endif

#if 0
// @return namespace prefix if `str` is prefixed or not-prefixed `target`
[[deprecated]] boost::optional<std::string> xmlnsXXX(const std::string& target,
                                                     const std::string& str) {
  if (str == target) {  // not prefixed
    return std::string{""};
  }
  std::smatch m;
  if (std::regex_match(
          str, m,
          std::regex{(boost::format{R"((?:(\w*):)?%1%)"} % target).str()})) {
    return {m.str(1)};
  }
  return {};
}

[[deprecated]] std::unordered_map<std::string, std::string> xmlattr(
    const boost::property_tree::ptree& pt) {
  std::unordered_map<std::string, std::string> r;
  for (const auto& it : pt) {
    if (it.first == "<xmlattr>") {
      for (const auto& i : it.second) {
        r[i.first] = i.second.data();
      }
    }
  }
  return r;
}
#endif


  // http://www.xml-cml.org/convention/molecular#atom-x2
  // > An atom MAY have an x2 attribute, the value of which is used for
  // > displaying the object in 2 dimensions. This is unrelated to the 3-D
  // > coordinates for the object.
  const auto x2 = atom.get_optional<double>(
      "<xmlattr>.x2", LexicalTranslator<double>{xpath_to_atom + "/@x2"});
  const auto y2 = atom.get_optional<double>(
      "<xmlattr>.y2", LexicalTranslator<double>{xpath_to_atom + "/@x2"});
 else if (x2 || y2) {
    if (!(x2 && y2)) {
      // http://www.xml-cml.org/convention/molecular#atom-x2
      // > If a x2 attribute is present there MUST also be a y2 attribute.
      auto msg =
          boost::format{"%1% does not have both of x2 and y2 attributes"} %
          xpath_to_atom;
      throw RDKit::FileParseException{msg.str()};
    }
    RDGeom::Point3D r{*x2, *y2, 0.0};
    BOOST_LOG(rdDebugLog) << xpath_to_atom << ' ' << r << std::endl;
    conformer->setAtomPos(idx, r);
  }

  if (spinMultiplicity && *spinMultiplicity != overall_spinMultiplicity) {
    auto msg =
        boost::format{
            "%1% (= %2%) is not equal to "
            "%3% = %4% (sum of %5%) - %6% (number of atoms)"} %
        xpath("@spinMultiplicity") % *spinMultiplicity %
        overall_spinMultiplicity % (overall_spinMultiplicity + num_atoms) %
        "../atomArray/atom/@spinMultiplicity" % num_atoms;
    throw RDKit::FileParseException{msg.str()};
  }
