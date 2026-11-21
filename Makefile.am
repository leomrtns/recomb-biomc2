SUBDIRS = src
dist_bin_SCRIPTS = scripts/biomc2.plot
EXTRA_DIST = README.md doc example
EXTRA_DIST += tests/regression.py conda/meta.yaml conda/build.sh
check-local:
	python3 $(top_srcdir)/tests/regression.py
