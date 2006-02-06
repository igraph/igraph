
EXTRA_DIST = $(TESTSUITE_AT) testsuite

MAINTAINERCLEANFILES = Makefile.in $(TESTSUITE) package.m4 atconfig

$(srcdir)/package.m4: $(top_srcdir)/configure.in
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
	visitors.at layout.at motifs.at topology.at
TESTSUITE = $(srcdir)/testsuite

AUTOTEST = $(AUTOM4TE) --language=autotest
$(srcdir)/$(TESTSUITE): package.m4 $(TESTSUITE_AT)
	$(AUTOTEST) -I $(srcdir) testsuite.at -o $@.tmp
	mv $@.tmp $@

clean-local: $(TESTSUITE)
	$(SHELL) $(TESTSUITE) --clean

check-local: atconfig atlocal $(TESTSUITE)
	$(SHELL) $(TESTSUITE)

# Run the test suite on the *installed* tree.
installcheck-local:
	$(SHELL) $(TESTSUITE) AUTOTEST_PATH=$(exec_prefix)/bin
