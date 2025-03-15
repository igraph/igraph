CXXFLAGS += -Wall -O3 -std=c++14
LDFLAGS += 

# Set INFOMAP_DIR to your Infomap directory
INFOMAP_DIR = ../../..
INFOMAP_LIB = $(INFOMAP_DIR)/lib/libInfomap.a

.PHONY: clean distclean

example: example.cpp $(INFOMAP_LIB) Makefile
	$(CXX) $(CXXFLAGS) -DNS_INFOMAP $< -o $@ -I$(INFOMAP_DIR)/include -L$(INFOMAP_DIR)/lib -lInfomap

mem-example: mem-example.cpp $(INFOMAP_LIB) Makefile
	$(CXX) $(CXXFLAGS) -DNS_INFOMAP $< -o $@ -I$(INFOMAP_DIR)/include -L$(INFOMAP_DIR)/lib -lInfomap

multi-example: multi-example.cpp $(INFOMAP_LIB) Makefile
	$(CXX) $(CXXFLAGS) -DNS_INFOMAP $< -o $@ -I$(INFOMAP_DIR)/include -L$(INFOMAP_DIR)/lib -lInfomap

file-example: file-example.cpp $(INFOMAP_LIB) Makefile
	$(CXX) $(CXXFLAGS) -DNS_INFOMAP $< -o $@ -I$(INFOMAP_DIR)/include -L$(INFOMAP_DIR)/lib -lInfomap

meta-example: meta-example.cpp $(INFOMAP_LIB) Makefile
	$(CXX) $(CXXFLAGS) -DNS_INFOMAP $< -o $@ -I$(INFOMAP_DIR)/include -L$(INFOMAP_DIR)/lib -lInfomap

$(INFOMAP_LIB):
	$(MAKE) -C $(INFOMAP_DIR) lib

clean:
	$(RM) example

distclean:
	$(MAKE) -C $(INFOMAP_DIR) clean