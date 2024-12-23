CXXFLAGS = -Wall -O3

# Set INFOMAP_DIR to your Infomap directory
INFOMAP_DIR = ../..
R_BUILD_DIR = $(INFOMAP_DIR)/build/R
MODULE_FILES = $(R_BUILD_DIR)/infomap.R $(R_BUILD_DIR)/infomap_wrap.so

.PHONY: clean distclean

infomap/infomap.R: $(R_BUILD_DIR)/infomap.R Makefile
	@echo Copying Infomap library and the R wrapper to a local folder...
	@mkdir -p infomap
	@cp -av $(R_BUILD_DIR)/infomap.R infomap
	@cp -av $(R_BUILD_DIR)/infomap_wrap.so infomap/infomap.so

$(R_BUILD_DIR)/infomap.R:
	cd $(INFOMAP_DIR) && $(MAKE) R

clean:
	$(RM) -r infomap

distclean:
	cd $(INFOMAP_DIR) && $(MAKE) clean