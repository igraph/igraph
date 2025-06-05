CXXFLAGS += -Wall -Wextra -pedantic -Wnon-virtual-dtor -std=c++14
LDFLAGS +=
CXX_CLANG := $(shell $(CXX) --version 2>/dev/null | grep clang)
BREW := $(shell which brew 2>/dev/null)

ifneq ($(BREW),)
	CXXFLAGS += -I$(shell brew --prefix)/include
	LDFLAGS += -L$(shell brew --prefix)/lib
endif

ifeq "$(findstring debug, $(MAKECMDGOALS))" "debug"
	CXXFLAGS += -O0 -g
else
	ifeq "$(CXX_CLANG)" ""
		CXXFLAGS += -O4
		ifneq "$(findstring noomp, $(MAKECMDGOALS))" "noomp"
			CXXFLAGS += -fopenmp
			LDFLAGS += -fopenmp
		endif
	else
		CXXFLAGS += -Wshadow -O3
		ifneq "$(findstring noomp, $(MAKECMDGOALS))" "noomp"
			CXXFLAGS += -Xpreprocessor -fopenmp
			LDFLAGS += -lomp
		endif
	endif
endif

##################################################
# General file dependencies
##################################################

HEADERS := $(shell find src -name "*.h")
SOURCES := $(shell find src -name "*.cpp")
OBJECTS := $(SOURCES:src/%.cpp=build/Infomap/%.o)

##################################################
# Stand-alone C++ targets
##################################################

.PHONY: all noomp test debug format

all: Infomap
	@true

Infomap: $(OBJECTS)
	@echo "Linking object files to target $@..."
	$(CXX) $(LDFLAGS) -o $@ $^
	@echo "-- Link finished --"

## Generic compilation rule for object files from cpp files
build/Infomap/%.o : src/%.cpp $(HEADERS) Makefile
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

noomp: Infomap
	@true

debug: Infomap
	@true

format: js-format py-format
	clang-format -i $(HEADERS) $(SOURCES)

##################################################
# JavaScript through Emscripten
##################################################

.PHONY: js-worker js-clean js-test

WORKER_FILENAME := infomap.worker.js
PRE_WORKER_MODULE := interfaces/js/pre-worker-module.js

js-worker: build/js/$(WORKER_FILENAME) Infomap
	@echo "Built $^"
	@mkdir -p interfaces/js/src/worker
	cp build/js/* interfaces/js/src/worker/
	npm run build
	cp interfaces/js/README.md .

js-test:
	$(RM) -r package mapequation-infomap-*.tgz
	npm pack
	tar -xzvf mapequation-infomap-*.tgz
	cp package/index.js examples/js
	sed -i.'backup' -e 's/src=".*"/src="index.js"/' \
		examples/js/infomap-worker.html
	open examples/js/infomap-worker.html
	sleep 5
	$(RM) -r package
	$(RM) -r mapequation-infomap-*.tgz
	$(RM) -r examples/js/index.js
	$(RM) -r examples/js/infomap-worker.html
	mv examples/js/infomap-worker.html{.backup,}

build/js/infomap.worker.js: $(SOURCES) $(PRE_WORKER_MODULE)
	@echo "Compiling Infomap to run in a worker in the browser..."
	@mkdir -p $(dir $@)
	em++ -std=c++14 -O3 -s WASM=0 -s ALLOW_MEMORY_GROWTH=1 -s DISABLE_EXCEPTION_CATCHING=0 -s ENVIRONMENT=worker --pre-js $(PRE_WORKER_MODULE) -o build/js/$(WORKER_FILENAME) $(SOURCES)

js-clean:
	$(RM) -r build/js interfaces/js/src/worker index.js *.d.ts README.md

js-format:
	prettier --write interfaces/js

##################################################
# Static C++ library
##################################################

# Use separate object files to compile with definitions
# NS_INFOMAP: Wrap code in namespace infomap
# AS_LIB: Skip main function
LIB_DIR = build/lib
LIB_HEADERS := $(HEADERS:src/%.h=include/%.h)
LIB_OBJECTS := $(SOURCES:src/%.cpp=build/lib/%.o)

.PHONY: lib

lib: lib/libInfomap.a $(LIB_HEADERS)
	@echo "Wrote static library to lib/ and headers to include/"

lib/libInfomap.a: $(LIB_OBJECTS) Makefile
	@echo "Creating static library..."
	@mkdir -p lib
	ar rcs $@ $^

# Rule for $(LIB_HEADERS)
include/%.h: src/%.h
	@mkdir -p $(dir $@)
	@cp -a $^ $@

# Rule for $(LIB_OBJECTS)
build/lib/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -DNS_INFOMAP -DAS_LIB -c $< -o $@


##################################################
# General SWIG helpers
##################################################

SWIG_FILES := $(shell find interfaces/swig -name "*.i")


##################################################
# Python module
##################################################

PY_BUILD_DIR = build/py
PY_ONLY_HEADERS := $(HEADERS:%.h=$(PY_BUILD_DIR)/headers/%.h)
PY_HEADERS := $(HEADERS:src/%.h=$(PY_BUILD_DIR)/src/%.h)
PY_SOURCES := $(SOURCES:src/%.cpp=$(PY_BUILD_DIR)/src/%.cpp)

.PHONY: python py-swig py-build

# Use python distutils to compile the module
python: py-swig py-build Makefile
	@true

py-build: Makefile
	@cp -a interfaces/python/setup.py $(PY_BUILD_DIR)/
	@touch $(PY_BUILD_DIR)/__init__.py
	@python utils/create-python-package-meta.py $(PY_BUILD_DIR)/package_meta.py
	@cat $(PY_BUILD_DIR)/package_meta.py $(PY_BUILD_DIR)/infomap.py > $(PY_BUILD_DIR)/temp.py
	@mv $(PY_BUILD_DIR)/temp.py $(PY_BUILD_DIR)/infomap.py
	@autopep8 --jobs 8 --aggressive --aggressive -i $(PY_BUILD_DIR)/infomap.py
	@cp -a interfaces/python/MANIFEST.in $(PY_BUILD_DIR)/
	@cp -a README.rst $(PY_BUILD_DIR)/
	@cp -a LICENSE_AGPLv3.txt $(PY_BUILD_DIR)/LICENSE
	@cd $(PY_BUILD_DIR) && CC=$(CXX) python3 setup.py build_ext --inplace

# Generate wrapper files from source and interface files
py-swig: $(PY_HEADERS) $(PY_SOURCES) $(PY_ONLY_HEADERS) interfaces/python/infomap.py
	@mkdir -p $(PY_BUILD_DIR)
	@cp -a $(SWIG_FILES) $(PY_BUILD_DIR)/
	swig -c++ -python -outdir $(PY_BUILD_DIR) -o $(PY_BUILD_DIR)/infomap_wrap.cpp $(PY_BUILD_DIR)/Infomap.i

# Rule for $(PY_HEADERS) and $(PY_SOURCES)
$(PY_BUILD_DIR)/src/%: src/%
	@mkdir -p $(dir $@)
	@cp -a $^ $@

$(PY_BUILD_DIR)/headers/%: %
	@mkdir -p $(dir $@)
	@cp -a $^ $@

.PHONY: py-doc py-local-install
SPHINX_SOURCE_DIR = interfaces/python/source
SPHINX_TARGET_DIR = docs

py-test:
	@cp -r examples/networks/*.net $(PY_BUILD_DIR)
	python3 -m flake8 --count --show-source --statistics --ignore E501,F811,W503 $(PY_BUILD_DIR)/infomap.py
	cd $(PY_BUILD_DIR) && python3 -m doctest infomap.py
	cd examples/python && for f in *.py; do python3 "$$f" > /dev/null || exit 1; done

py-local-install:
	# Run this to get 'import infomap' to always import the latest
	# locally built version, so no need to run this multiple times.
	pip install -e $(PY_BUILD_DIR)

py-doc: py-local-install
	# Uses docstrings from the infomap available with 'import infomap'.
	# Run py-local-install if you don't have pip installed it with -e
	# and don't have the latest version installed
	@mkdir -p $(SPHINX_TARGET_DIR)
	@cp -a README.rst ${SPHINX_SOURCE_DIR}/index.rst
	sphinx-build -b html $(SPHINX_SOURCE_DIR) $(SPHINX_TARGET_DIR)
	@rm -r ${SPHINX_SOURCE_DIR}/index.rst
	npx prettier --write docs/searchindex.js

.PHONY: pypitest-publish pypi-publish py-clean
PYPI_DIR = $(PY_BUILD_DIR)
PYPI_SDIST = $(shell find $(PYPI_DIR) -name "*.tar.gz" 2>/dev/null)

py-clean:
	$(RM) -r $(PY_BUILD_DIR)/dist

py-format:
	python3 -m isort interfaces/python examples/python
	python3 -m black interfaces/python examples/python

pypi-dist:
	cd $(PY_BUILD_DIR) && python setup.py sdist bdist_wheel

# pip -vvv --no-cache-dir install --upgrade -I --index-url https://test.pypi.org/simple/ infomap
# pip install -e build/py/pypi/infomap/
pypitest-publish:
	@[ "${PYPI_SDIST}" ] && echo "Publish dist..." || ( echo "dist files not built"; exit 1 )
	cd $(PYPI_DIR) && python -m twine upload -r testpypi --verbose dist/*

pypi-publish:
	@[ "${PYPI_SDIST}" ] && echo "Publish dist..." || ( echo "dist files not built"; exit 1 )
	cd $(PYPI_DIR) && python -m twine upload --skip-existing --verbose dist/*


##################################################
# R module
##################################################

R_BUILD_DIR = build/R
R_HEADERS := $(HEADERS:src/%.h=$(R_BUILD_DIR)/src/%.h)
R_SOURCES := $(SOURCES:src/%.cpp=$(R_BUILD_DIR)/src/%.cpp)

.PHONY: R R-build

# Use R to compile the module
R: R-build Makefile
	@mkdir -p $(R_BUILD_DIR)/Infomap
	@cp -a examples/R/load-infomap.R $(R_BUILD_DIR)/Infomap/
	cd $(R_BUILD_DIR) && CXX="$(CXX)" PKG_CPPFLAGS="$(CXXFLAGS) -DAS_LIB" PKG_LIBS="$(LDFLAGS)" R CMD SHLIB infomap_wrap.cpp $(SOURCES)
	@cp -a $(R_BUILD_DIR)/infomap.R $(R_BUILD_DIR)/Infomap/
	@cp -a $(R_BUILD_DIR)/infomap_wrap.so $(R_BUILD_DIR)/Infomap/infomap.so
	@cp -a examples/R/example-minimal.R $(R_BUILD_DIR)/Infomap/
	@true

# Generate wrapper files from source and interface files
R-build: Makefile $(R_HEADERS) $(R_SOURCES)
	@mkdir -p $(R_BUILD_DIR)
	@cp -a $(SWIG_FILES) $(R_BUILD_DIR)/
	swig -c++ -r -outdir $(R_BUILD_DIR) -o $(R_BUILD_DIR)/infomap_wrap.cpp $(R_BUILD_DIR)/Infomap.i

# Rule for $(R_HEADERS) and $(R_SOURCES)
$(R_BUILD_DIR)/src/%: src/%
	@mkdir -p $(dir $@)
	@cp -a $^ $@


##################################################
# Docker
##################################################

TAG_NAME = mapequation/infomap

docker-build: Makefile docker/infomap.Dockerfile
	docker-compose build infomap

docker-run: Makefile
	docker-compose run --rm infomap

docker-build-notebook: Makefile docker/notebook.Dockerfile
	docker-compose build notebook

docker-run-notebook: Makefile
	docker-compose up notebook

# rstudio
docker-build-rstudio: Makefile docker/rstudio.Dockerfile
	docker build \
	-f docker/rstudio.Dockerfile \
	-t $(TAG_NAME):rstudio .

docker-run-rstudio: Makefile
	docker run --rm \
	$(TAG_NAME):rstudio

# ubuntu test python
docker-build-ubuntu-test-python: Makefile
	docker build -f docker/ubuntu.Dockerfile -t infomap:python-test .

docker-run-ubuntu-test-python: Makefile
	docker run --rm infomap:python-test

# R with RStudio
docker-build-r: Makefile
	docker build -f docker/rstudio.Dockerfile -t infomap:r .

docker-run-r: Makefile
	docker run --rm -p 8787:8787 -e PASSWORD=InfomapR infomap:r

##################################################
# Clean
##################################################

clean: js-clean
	$(RM) -r Infomap build lib include
