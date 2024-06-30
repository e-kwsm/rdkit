//
//  Copyright (C) 2022-2023 Tad Hurst, Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

//  This software is based on the Chemaxon documentation for the MRV format:P
// https://docs.chemaxon.com/display/docs/marvin-documents-mrv.md
// and this implmentation is tested against the parsing and generation in the
// Marvin JS sketcher: https://marvinjs-demo.chemaxon.com/latest/demo.html

#ifndef RD_MARVINDEFS_H
#define RD_MARVINDEFS_H

#include <GraphMol/RDKitBase.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/xml_parser.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <float.h>  // Needed for DBL_MAX on Clang

using boost::property_tree::ptree;

namespace RDKit {

const std::vector<std::string> sruSgroupConnectChoices{"hh", "ht", "eu"};
const std::vector<std::string> marvinBondOrders{"1", "2", "3", "A"};
const std::vector<std::string> marvinQueryBondsTypes{"SD", "SA", "DA", "Any"};
const std::vector<std::string> marvinConventionTypes{"cxn:coord",
                                                     "cxn:hydrogen"};
const std::vector<std::string> marvinStereoDictRefTypes{"cml:W", "cml:H"};
const std::vector<std::string> marvinStereoConventionTypes{"1", "3", "4", "6"};

const std::vector<std::string> marvinRadicalVals{
    "monovalent", "divalent",   "divalent1",  "divalent3",
    "trivalent",  "trivalent2", "trivalent4", "4"};
const std::map<std::string, int> marvinRadicalToRadicalElectrons{
    {"monovalent", 1}, {"divalent", 2},   {"divalent1", 2},  {"divalent3", 2},
    {"trivalent", 3},  {"trivalent2", 3}, {"trivalent4", 3}, {"4", 4}};

const std::map<int, std::string> radicalElectronsToMarvinRadical{
    {1, "monovalent"}, {2, "divalent"}, {3, "trivalent4"}, {4, "4"}};

enum IsSgroupInAtomSetResult {
  SgroupInAtomSet,
  SgroupNotInAtomSet,
  SgroupBothInAndNotInAtomSet,
};

class MarvinWriterException : public std::runtime_error {
 public:
  explicit MarvinWriterException(std::string message)
      : std::runtime_error(message) {};
};

class MarvinArrow {
 public:
  std::string type;
  double x1;
  double y1;
  double x2;
  double y2;

  [[nodiscard]] std::string toString() const;
  [[nodiscard]] ptree toPtree() const;
};

class MarvinPlus {
 public:
  std::string id;
  double x1;
  double y1;
  double x2;
  double y2;

  [[nodiscard]] std::string toString() const;
  [[nodiscard]] ptree toPtree() const;
};

class MarvinCondition {
 public:
  std::string id;
  std::string text;
  double x;
  double y;
  double fontScale = 0.0;

  std::string halign;
  std::string valign;

  [[nodiscard]] std::string toString() const;
  [[nodiscard]] ptree toPtree() const;
};

class MarvinAttachmentPoint {
 public:
  // <attachmentPoint atom="a7" order="1" bond="b6"/>
  std::string atom;
  std::string bond;
  std::string order;

  [[nodiscard]] std::string toString() const;
  [[nodiscard]] ptree toPtree() const;
};

class MarvinAtom {
 public:
  std::string id;
  std::string elementType;
  double x2;
  double y2;
  double x3;
  double y3;
  double z3;

  int formalCharge;
  std::string radical;
  int isotope;
  int mrvValence;
  int hydrogenCount;
  std::string mrvAlias;
  std::string mrvStereoGroup;
  int mrvMap;
  std::string sgroupRef;
  bool sGroupRefIsSuperatom;  // if set, we will not change the sgroupRef - the
                              // superatom really needs it
  std::string sgroupAttachmentPoint;
  int rgroupRef;

  MarvinAtom();
  MarvinAtom(const MarvinAtom &atomToCopy, std::string newId);

  bool operator==(const MarvinAtom &rhs) const;

  bool operator==(const MarvinAtom *rhs) const;

  [[nodiscard]] bool isElement() const;

  [[nodiscard]] std::string toString() const;
  [[nodiscard]] ptree toPtree(unsigned int coordinatePrecision = 6) const;
};

class MarvinBondStereo {
 public:
  std::string value;
  std::string convention;
  std::string conventionValue;
  std::string dictRef;

  [[nodiscard]] std::string toString() const;
  [[nodiscard]] ptree toPtree() const;
};

class MarvinBond {
 public:
  std::string id;
  std::string atomRefs2[2];
  std::string order;
  MarvinBondStereo bondStereo;
  std::string queryType;
  std::string convention;

  MarvinBond() {}

  MarvinBond(const MarvinBond &bondToCopy, std::string newId,
             std::string atomRef1, std::string atomRef2);

  [[nodiscard]] bool isEqual(const MarvinAtom &other) const;

  bool operator==(const MarvinAtom &rhs) const;

  [[nodiscard]] const std::string getBondType() const;

  [[nodiscard]] std::string toString() const;
  [[nodiscard]] ptree toPtree() const;
};

class MarvinRectangle {
 protected:
  RDGeom::Point3D center;
  bool centerIsStale = true;

 public:
  RDGeom::Point3D upperLeft;
  RDGeom::Point3D lowerRight;

  MarvinRectangle(double left, double right, double top, double bottom);
  MarvinRectangle(const RDGeom::Point3D &upperLeftInit,
                  const RDGeom::Point3D &lowerRightInit);
  MarvinRectangle(const std::vector<MarvinAtom *> &atoms);
  MarvinRectangle(const std::vector<MarvinRectangle> &rects);

  void extend(const MarvinRectangle &otherRectangle);

  RDGeom::Point3D &getCenter();

  [[nodiscard]] bool overlapsVertically(
      const MarvinRectangle &otherRectangle) const;

  [[nodiscard]] bool overlapsVHorizontally(
      const MarvinRectangle &otherRectangle) const;

  static bool compareRectanglesByX(MarvinRectangle &r1, MarvinRectangle &r2);

  static bool compareRectanglesByYReverse(MarvinRectangle &r1,
                                          MarvinRectangle &r2);
};

template <typename T>
typename std::vector<std::unique_ptr<T>>::iterator findUniquePtr(
    std::vector<std ::unique_ptr<T>> &vector, T *itemToFind) {
  auto foundUniqIter = find_if(vector.begin(), vector.end(),
                               [itemToFind](std::unique_ptr<T> &uniquePtr) {
                                 return uniquePtr.get() == itemToFind;
                               });

  if (foundUniqIter == vector.end()) {
    throw FileParseException("Unexpected error - item to find not found");
  }

  return foundUniqIter;
}

template <typename T>
void eraseUniquePtr(std::vector<std ::unique_ptr<T>> &vector, T *itemToErase) {
  auto removeUniqIter = findUniquePtr<T>(vector, itemToErase);

  if (removeUniqIter == vector.end()) {
    throw FileParseException("Unexpected error - item to remove not found");
  }

  vector.erase(removeUniqIter);
}

class MarvinMolBase {
 public:
  std::string molID;
  std::string id;  // used in all sGroups
  unsigned int coordinatePrecision = 6;
  std::vector<MarvinAtom *> atoms;  // owned by parent MarvinMol
  std::vector<MarvinBond *> bonds;  // owned by parent MarvinMol
  std::vector<std::unique_ptr<MarvinMolBase>> sgroups;
  MarvinMolBase *parent;

  [[nodiscard]] virtual std::string role() const = 0;
  [[nodiscard]] virtual bool hasAtomBondBlocks() const = 0;
  [[nodiscard]] virtual std::string toString() const = 0;
  [[nodiscard]] virtual ptree toPtree() const;
  void addSgroupsToPtree(ptree &pt) const;

  [[nodiscard]] virtual MarvinMolBase *copyMol(
      const std::string &idAppend) const = 0;
  virtual void pushOwnedAtom(MarvinAtom *atom);
  virtual void pushOwnedBond(MarvinBond *bond);

  virtual void removeOwnedAtom(MarvinAtom *atom);
  virtual void removeOwnedBond(MarvinBond *bond);

  void setPrecision(unsigned int precision);

  [[nodiscard]] int getExplicitValence(const MarvinAtom &marvinAtom) const;

  MarvinMolBase() {}

  virtual ~MarvinMolBase();

  [[nodiscard]] int getAtomIndex(std::string id) const;
  [[nodiscard]] int getBondIndex(std::string id) const;

  [[nodiscard]] const std::vector<std::string> getBondList() const;
  [[nodiscard]] const std::vector<std::string> getAtomList() const;
  bool AnyOverLappingAtoms(const MarvinMolBase *otherMol) const;

  void cleanUpNumbering(
      int &molCount  // this is the starting mol count, and receives the ending
                     // mol count - THis is used when
                     // MarvinMol->convertToSuperAtoms is called multiple times
                     // from a RXN
      ,
      int &atomCount  // starting and ending atom count
      ,
      int &bondCount  // starting and ending bond count
      ,
      int &sgCount  // starting and ending sq count
      ,
      std::map<std::string, std::string>
          &sgMap  // map from old sg number to new sg number
      ,
      std::map<std::string, std::string>
          &atomMap  // map from old atom number to new atom number
      ,
      std::map<std::string, std::string>
          &bondMap  // map from old bond number to new bond number
  );

  // the following is virtual because some derived classes need to do more than
  // just call the base class.  Currently, only MarvinSuperatomSgroup does this
 public:
  virtual void cleanUpNumberingMolsAtomsBonds(
      int &molCount,  // this is the starting mol count, and receives the ending
                      // mol count - THis is used when
                      // MarvinMol->convertToSuperAtoms is called multiple
                      // times from a RXN
      int &atomCount,  // starting and ending atom count
      int &bondCount,  // starting and ending bond count
      std::map<std::string, std::string> &sgMap,
      std::map<std::string, std::string> &atomMap,
      std::map<std::string, std::string> &bondMap);

  void cleanUpSgNumbering(int &sgCount,
                          std::map<std::string, std::string> &sgMap);

  // the following is virtual because some derived classes need to do more than
  // just call the base class.  Currently, only MarvinSuperatomSgroup does this

  [[nodiscard]] virtual IsSgroupInAtomSetResult isSgroupInSetOfAtoms(
      const std::vector<MarvinAtom *> &setOfAtoms) const;

 public:
  static bool atomRefInAtoms(MarvinAtom *a, std::string b);
  static bool bondRefInBonds(MarvinBond *a, std::string b);
  static bool molIDInSgroups(std::string a, std::string b);
  MarvinAtom *findAtomByRef(std::string atomId);
  MarvinBond *findBondByRef(std::string atomId);

  void prepSgroupsForRDKit();
  void processSgroupsFromRDKit();

  [[nodiscard]] virtual bool isPassiveRoleForExpansion() const;
  [[nodiscard]] virtual bool isPassiveRoleForContraction() const;
  virtual void processSpecialSgroups();
  virtual void parseMoleculeSpecific(RDKit::RWMol *mol,
                                     std::unique_ptr<SubstanceGroup> &sgroup,
                                     int sequenceId);

  [[nodiscard]] bool has2dCoords() const;
  [[nodiscard]] bool has3dCoords() const;
  [[nodiscard]] bool hasAny3dCoords() const;
  [[nodiscard]] bool hasAny2dCoords() const;
  [[nodiscard]] bool hasCoords() const;
  void removeCoords();

  void parseAtomsAndBonds(ptree &molTree);
};

class MarvinSruCoModSgroup : public MarvinMolBase {
 private:
  std::string roleName;  // could be MarvinSruSgroup, MarvinCopolymerSgroup or
                         // MarvinModificationSgroup
 public:
  [[nodiscard]] MarvinMolBase *copyMol(
      const std::string &idAppend) const override;
  MarvinSruCoModSgroup(std::string type, MarvinMolBase *parent);
  MarvinSruCoModSgroup(MarvinMolBase *parent, std::string role, ptree &molTree);

  std::string title;
  std::string connect;
  std::string correspondence;

  [[nodiscard]] std::string toString() const override;
  [[nodiscard]] ptree toPtree() const override;

  [[nodiscard]] std::string role() const override;
  [[nodiscard]] bool hasAtomBondBlocks() const override;
  void parseMoleculeSpecific(RDKit::RWMol *mol,
                             std::unique_ptr<SubstanceGroup> &sgroup,
                             int sequenceId) override;
};

class MarvinDataSgroup : public MarvinMolBase {
 public:
  MarvinDataSgroup(MarvinMolBase *parent);
  MarvinDataSgroup(MarvinMolBase *parent, ptree &molTree);

  [[nodiscard]] MarvinMolBase *copyMol(
      const std::string &idAppend) const override;

  std::string context;
  std::string fieldName;
  std::string placement;
  std::string unitsDisplayed;
  std::string queryType;
  std::string queryOp;
  std::string fieldData;
  std::string units;
  double x;
  double y;

  [[nodiscard]] std::string toString() const override;
  [[nodiscard]] ptree toPtree() const override;

  [[nodiscard]] std::string role() const override;
  [[nodiscard]] bool hasAtomBondBlocks() const override;
  void parseMoleculeSpecific(RDKit::RWMol *mol,
                             std::unique_ptr<SubstanceGroup> &sgroup,
                             int sequenceId) override;
};

class MarvinSuperatomSgroupExpanded : public MarvinMolBase {
 public:
  std::string title;

  MarvinSuperatomSgroupExpanded(MarvinMolBase *parent);
  MarvinSuperatomSgroupExpanded(MarvinMolBase *parent, ptree &molTree);

  [[nodiscard]] MarvinMolBase *copyMol(
      const std::string &idAppend) const override;

  ~MarvinSuperatomSgroupExpanded() override;

  MarvinMolBase *convertToOneSuperAtom();

  [[nodiscard]] std::string toString() const override;
  [[nodiscard]] ptree toPtree() const override;

  [[nodiscard]] std::string role() const override;
  [[nodiscard]] bool hasAtomBondBlocks() const override;
  [[nodiscard]] bool isPassiveRoleForContraction() const override;

  void parseMoleculeSpecific(RDKit::RWMol *mol,
                             std::unique_ptr<SubstanceGroup> &sgroup,
                             int sequenceId) override;
};

class MarvinMultipleSgroup : public MarvinMolBase {
 public:
  MarvinMultipleSgroup(MarvinMolBase *parent);
  MarvinMultipleSgroup(MarvinMolBase *parent, ptree &molTree);

  [[nodiscard]] MarvinMolBase *copyMol(
      const std::string &idAppend) const override;

  std::string title;
  bool isExpanded = false;
  std::vector<MarvinAtom *> parentAtoms;
  std::vector<MarvinBond *>
      bondsToAtomsNotInExpandedGroup;  // only when expanded

  void expandOneMultipleSgroup();
  void contractOneMultipleSgroup();
  int getMatchedOrphanBondIndex(std::string atomIdToCheck,
                                std::vector<MarvinBond *> &bondsToTry,
                                std::vector<MarvinBond *> &orphanedBonds) const;

  [[nodiscard]] std::string toString() const override;
  [[nodiscard]] ptree toPtree() const override;

  [[nodiscard]] std::string role() const override;
  [[nodiscard]] bool hasAtomBondBlocks() const override;
  [[nodiscard]] bool isPassiveRoleForExpansion() const override;
  [[nodiscard]] bool isPassiveRoleForContraction() const override;
  void processSpecialSgroups() override;

  void parseMoleculeSpecific(RDKit::RWMol *mol,
                             std::unique_ptr<SubstanceGroup> &sgroup,
                             int sequenceId) override;
};

class MarvinMulticenterSgroup : public MarvinMolBase {
  // <molecule molID="m2" id="sg1" role="MulticenterSgroup" atomRefs="a2 a6 a5
  // a4 a3" center="a18"/>sgroup->
 public:
  MarvinMulticenterSgroup(MarvinMolBase *parent);
  MarvinMulticenterSgroup(MarvinMolBase *parent, ptree &molTree);

  [[nodiscard]] MarvinMolBase *copyMol(
      const std::string &idAppend) const override;

  void processOneMulticenterSgroup();

  [[nodiscard]] std::string toString() const override;
  [[nodiscard]] ptree toPtree() const override;

  MarvinAtom *center;
  [[nodiscard]] std::string role() const override;
  [[nodiscard]] bool hasAtomBondBlocks() const override;
  void processSpecialSgroups() override;
  void parseMoleculeSpecific(RDKit::RWMol *mol,
                             std::unique_ptr<SubstanceGroup> &sgroup,
                             int sequenceId) override;
};

class MarvinGenericSgroup : public MarvinMolBase {
  // <molecule molID="m2" id="sg1" role="GenericSgroup" atomRefs="a1 a2 a3 a4 a5
  // a6 a7 a8 a9 a13 a10 a11 a12" charge="onAtoms"/></molecule>
 public:
  MarvinGenericSgroup(MarvinMolBase *parent);
  MarvinGenericSgroup(MarvinMolBase *parent, ptree &molTree);

  [[nodiscard]] MarvinMolBase *copyMol(
      const std::string &idAppend) const override;

  std::string charge;  // onAtoms or onBrackets
  [[nodiscard]] std::string toString() const override;
  [[nodiscard]] ptree toPtree() const override;

  [[nodiscard]] std::string role() const override;
  [[nodiscard]] bool hasAtomBondBlocks() const override;
  void parseMoleculeSpecific(RDKit::RWMol *mol,
                             std::unique_ptr<SubstanceGroup> &sgroup,
                             int sequenceId) override;
};

class MarvinMonomerSgroup : public MarvinMolBase {
  // <molecule id="sg1" role="MonomerSgroup" title="mon" charge="onAtoms"
  // molID="m2" atomRefs="a2 a1 a3 a4">
  // </molecule>
 public:
  MarvinMonomerSgroup(MarvinMolBase *parent);
  MarvinMonomerSgroup(MarvinMolBase *parent, ptree &molTree);

  [[nodiscard]] MarvinMolBase *copyMol(
      const std::string &idAppend) const override;

  std::string title;
  std::string charge;  // onAtoms or onBrackets
  [[nodiscard]] std::string toString() const override;
  [[nodiscard]] ptree toPtree() const override;

  [[nodiscard]] std::string role() const override;
  [[nodiscard]] bool hasAtomBondBlocks() const override;
  void parseMoleculeSpecific(RDKit::RWMol *mol,
                             std::unique_ptr<SubstanceGroup> &sgroup,
                             int sequenceId) override;
};

class MarvinSuperatomSgroup : public MarvinMolBase {
 public:
  std::string title;
  std::vector<std::unique_ptr<MarvinAttachmentPoint>> attachmentPoints;

  MarvinSuperatomSgroup(MarvinMolBase *parent);
  MarvinSuperatomSgroup(MarvinMolBase *parent, ptree &molTree);

  [[nodiscard]] MarvinMolBase *copyMol(
      const std::string &idAppend) const override;

  ~MarvinSuperatomSgroup() override;

  void convertFromOneSuperAtom();

  [[nodiscard]] std::string role() const override;
  [[nodiscard]] bool hasAtomBondBlocks() const override;
  [[nodiscard]] bool isPassiveRoleForExpansion() const override;

  [[nodiscard]] std::string toString() const override;
  [[nodiscard]] ptree toPtree() const override;

  void cleanUpNumberingMolsAtomsBonds(
      int &molCount,  // this is the starting mol count, and receives the ending
                      // mol count - THis is used when
                      // MarvinMol->convertToSuperAtoms is called multiple
                      // times from a RXN
      int &atomCount,  // starting and ending atom count
      int &bondCount,  // starting and ending bond count
      std::map<std::string, std::string> &sgMap,
      std::map<std::string, std::string> &atomMap,
      std::map<std::string, std::string> &bondMap) override;

  [[nodiscard]] IsSgroupInAtomSetResult isSgroupInSetOfAtoms(
      const std::vector<MarvinAtom *> &setOfAtoms) const override;

  void processSpecialSgroups() override;
};

class MarvinMol : public MarvinMolBase {
 public:
  MarvinMol();
  MarvinMol(ptree &molTree);

  [[nodiscard]] MarvinMolBase *copyMol(
      const std::string &idAppend) const override;

  ~MarvinMol() override;

  std::vector<std::unique_ptr<MarvinAtom>> ownedAtoms;
  std::vector<std::unique_ptr<MarvinBond>> ownedBonds;

  void pushOwnedAtom(MarvinAtom *atom) override;
  void pushOwnedBond(MarvinBond *bond) override;

  void removeOwnedAtom(MarvinAtom *atom) override;
  void removeOwnedBond(MarvinBond *bond) override;

  [[nodiscard]] std::string role() const override;
  [[nodiscard]] bool hasAtomBondBlocks() const override;
  [[nodiscard]] bool isPassiveRoleForContraction() const override;

  [[nodiscard]] std::string toString() const override;
  [[nodiscard]] ptree toPtree() const override;

  std::string generateMolString();
  [[nodiscard]] ptree toMolPtree() const;
};

class MarvinReaction {
 public:
  std::vector<std::unique_ptr<MarvinMol>> reactants;
  std::vector<std::unique_ptr<MarvinMol>> agents;
  std::vector<std::unique_ptr<MarvinMol>> products;

  MarvinArrow arrow;
  std::vector<std::unique_ptr<MarvinPlus>> pluses;
  std::vector<std::unique_ptr<MarvinCondition>> conditions;

  ~MarvinReaction();

  void prepSgroupsForRDKit();

  std::string toString();
  [[nodiscard]] ptree toPtree() const;
};

class MarvinStereoGroup {
 public:
  StereoGroupType groupType;  // one of ABS AND OR
  int groupNumber;
  std::vector<unsigned int> atoms;

  MarvinStereoGroup(StereoGroupType grouptypeInit, int groupNumberInit);
};

template <typename T>
bool getCleanNumber(std::string strToParse, T &outInt);
}  // namespace RDKit

#endif  // RD_MARVINDEFS_H
