ifneq ($(shell which clang++ 2> /dev/null),)
	CXX := clang++
else
	CXX := g++
endif

Boost_DIR := ${HOME}/opt/boost/latest

CXXFLAGS := -std=c++17 -Wall -Wextra -Wpedantic -g -fPIC
CPPFLAGS := \
	-isystem $(shell python3 -c 'import sysconfig;print(sysconfig.get_config_var("INCLUDEPY"))') \
	-isystem $(Boost_DIR)/include \
	-isystem ../Code \
	-isystem ../build/Code \
	-isystem ../External/catch/catch/single_include/catch2
LDFLAGS := \
	-L$(Boost_DIR)/lib -Wl,-rpath=$(Boost_DIR)/lib \
	-L../LLVM/lib -Wl,-rpath=../LLVM/lib
LDLIBS := -lRDKitFileParsers -lRDKitGraphMol -lRDKitRDGeometryLib -lRDKitRDGeneral \
	-lboost_python$(shell python3 -c 'import sysconfig;print(sysconfig.get_config_var("py_version_nodot"))') \
	$(shell python3 -c 'import sysconfig;print(sysconfig.get_config_var("BLDLIBRARY"))') \
	-lpthread -Wl,--as-needed

CMAKE_CXX_CLANG_TIDY := clang-tidy --fix --header-filter=! \
	--checks=bugprone-*,clang-*,misc-*,-misc-non-private-member-variables-in-classes,-misc-unused-parameters,modernize-*,-modernize-use-trailing-return-type,readbility-*

.PHONY: all
all: a.out file_parsers_catch.out rdmolfiles.so
a.out: main.o CMLReader.o
	$(CXX) -o$@ $^ $(LDFLAGS) $(LDLIBS)

file_parsers_catch.out: file_parsers_catch.o CMLReader.o
	$(CXX) -o$@ $^ $(LDFLAGS) $(LDLIBS)

rdmolfiles.so: rdmolfiles.o CMLReader.o
	$(CXX) -shared -fPIC -o$@ $^ $(LDFLAGS) $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o$@ -c $^

.PHONY: clang-format
clang-format:
	clang-format -i $(shell realpath CMLReader.cpp)

.PHONY: clang-tidy
clang-tidy: $(patsubst %.cpp,%.tidy,CMLReader.cpp file_parsers_catch.cpp main.cpp rdmolfiles.cpp)

%.tidy: %.cpp
	$(CMAKE_CXX_CLANG_TIDY) $(shell realpath $^) -- $(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $(shell realpath $^)
	touch $@

.PHONY: clean
clean:
	echo *.so *.o *.out *.tidy | xargs rm -f
