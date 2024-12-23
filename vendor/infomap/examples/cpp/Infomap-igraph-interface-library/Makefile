CXXFLAGS = -Wall -O3 -std=c++14

# Set INFOMAP_DIR to your Infomap directory
INFOMAP_DIR = ../../..
INFOMAP_LIB = $(INFOMAP_DIR)/lib/libInfomap.a
INC_INFOMAP = -I$(INFOMAP_DIR)/include -L$(INFOMAP_DIR)/lib -lInfomap

# igraph paths
IGRAPH_LIB_DIR = /usr/local/lib
IGRAPH_INCLUDE_DIR = /usr/local/include/igraph
INC_IGRAPH = -I$(IGRAPH_INCLUDE_DIR) -L$(IGRAPH_LIB_DIR) -ligraph

# Abbreviation
III = Infomap-igraph-interface

.PHONY: clean distclean

lib/lib$(III).a: build/$(III).o include/$(III).h
	@echo "Creating static $(III) library..."
	@mkdir -p lib
	ar rcs $@ $^

# Rule for lib headers
include/%.h: src/%.h
	@mkdir -p $(dir $@)
	@cp -a $^ $@

# Rule for lib objects
build/%.o: src/%.cpp src/%.h $(INFOMAP_LIB) Makefile
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -DNS_INFOMAP -c $< -o $@ -I$(IGRAPH_INCLUDE_DIR) -I$(INFOMAP_DIR)/include

# Make Infomap library
$(INFOMAP_LIB):
	cd $(INFOMAP_DIR) && $(MAKE) lib

clean:
	$(RM) -r build lib include

distclean:
	cd $(INFOMAP_DIR) && $(MAKE) clean