AUTOMAKE_OPTIONS = gnu

ACLOCAL_AMFLAGS= -I m4

SUBDIRS = src/common src man tests

src: src/common

.PHONY: test	
test: check
	cd tests && $(MAKE) $(AM_MAKEFLAGS) test-local

