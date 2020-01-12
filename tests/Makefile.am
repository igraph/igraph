
EXTRA_DIST = $(TESTSUITE_AT) $(top_builddir)/tests/testsuite

MAINTAINERCLEANFILES = Makefile.in $(TESTSUITE) package.m4 atconfig

$(srcdir)/package.m4: $(top_srcdir)/configure.ac
	{					\
	  echo '# Signature of the current package.'; \
	  echo 'm4_define([AT_PACKAGE_NAME],      [@PACKAGE_NAME@])'; \
	  echo 'm4_define([AT_PACKAGE_TARNAME],   [@PACKAGE_TARNAME@])'; \
	  echo 'm4_define([AT_PACKAGE_VERSION],   [@PACKAGE_VERSION@])'; \
	  echo 'm4_define([AT_PACKAGE_STRING],    [@PACKAGE_STRING@])'; \
	  echo 'm4_define([AT_PACKAGE_BUGREPORT], [@PACKAGE_BUGREPORT@])'; \
	} >$(srcdir)/package.m4

EXTRA_DIST += package.m4

TESTSUITE_AT = \
	testsuite.at \
	types.at basic.at structure_generators.at \
	structural_properties.at iterators.at components.at \
	visitors.at layout.at motifs.at topology.at foreign.at operators.at \
	other.at foreign.at conversion.at flow.at community.at eigen.at \
	cliques.at attributes.at arpack.at bipartite.at centralization.at \
	version.at separators.at hrg.at microscopic.at mt.at random.at scg.at \
	matching.at qsort.at coloring.at embedding.at

TESTSUITE = testsuite

AUTOTEST = $(AUTOM4TE) --language=autotest
$(TESTSUITE): $(srcdir)/package.m4 $(TESTSUITE_AT)
	$(AUTOTEST) -I $(top_srcdir)/tests $(top_srcdir)/tests/testsuite.at -o $@.tmp
	mv $@.tmp $@

clean-local: $(TESTSUITE)
	$(SHELL) $(TESTSUITE) --clean

check-local: atconfig atlocal $(TESTSUITE)
	if [ ! -f $(TESTSUITE) ]; then cp $(top_srcdir)/tests/testsuite .; fi
	$(SHELL) $(TESTSUITE)

# Run the test suite on the *installed* tree.
installcheck-local:
	$(SHELL) $(TESTSUITE) AUTOTEST_PATH=$(exec_prefix)/bin
