
rdkit_library(GeneralizedSubstruct 
              XQMol.cpp TextIO.cpp
              LINK_LIBRARIES MolEnumerator TautomerQuery SubstructMatch SmilesParse GraphMol Fingerprints)
target_compile_definitions(GeneralizedSubstruct PRIVATE RDKIT_GENERALIZEDSUBSTRUCT_BUILD)

rdkit_headers(XQMol.h DEST GraphMol/GeneralizedSubstruct)

if(RDK_USE_BOOST_SERIALIZATION)
rdkit_catch_test(generalizedSubstructCatch catch_tests.cpp LINK_LIBRARIES GenericGroups FileParsers GeneralizedSubstruct)
endif()

if (RDK_BUILD_PYTHON_WRAPPERS)
    add_subdirectory(Wrap)
endif ()
