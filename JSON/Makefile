CXX := clang++
CXXFLAGS := -Wall -Wextra -Wpedantic
CPPFLAGS := -I../Code -I../Code/GraphMol/FileParsers

LDFLAGS := -L../build/lib '-Wl,-rpath=$$ORIGIN/../build/lib'

LDLIBS := \
        -g \
        -lRDKitFileParsers \
        -lRDKitGraphMol \
        -lRDKitRDGeneral \
        -lRDKitMolHash \

.PHONY: all
all: zzz.cjson

a.out: a.o CJSONWriter.o
	$(CXX) -o $@ $^ $(LDFLAGS) $(LDLIBS)

CJSONWriter.o: CJson.h

zzz.json: a.out
	./$<

%.cjson: %.json
	bat $<
	python3 -m json.tool --indent 2 $< $@

.PHONY: clean
clean:
	rm -f CJSONWriter.o a.out
