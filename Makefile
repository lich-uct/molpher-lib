#
#  There exist several targets which are by default empty and which can be 
#  used for execution of your targets. These targets are usually executed 
#  before and after some main targets. They are: 
#
#     .build-pre:              called before 'build' target
#     .build-post:             called after 'build' target
#     .clean-pre:              called before 'clean' target
#     .clean-post:             called after 'clean' target
#     .clobber-pre:            called before 'clobber' target
#     .clobber-post:           called after 'clobber' target
#     .all-pre:                called before 'all' target
#     .all-post:               called after 'all' target
#     .help-pre:               called before 'help' target
#     .help-post:              called after 'help' target
#
#  Targets beginning with '.' are not intended to be called on their own.
#
#  Main targets can be executed directly, and they are:
#  
#     build                    build a specific configuration
#     clean                    remove built files from a configuration
#     clobber                  remove all built files
#     all                      build all configurations
#     help                     print help mesage
#  
#  Targets .build-impl, .clean-impl, .clobber-impl, .all-impl, and
#  .help-impl are implemented in nbproject/makefile-impl.mk.
#
#  Available make variables:
#
#     CND_BASEDIR                base directory for relative paths
#     CND_DISTDIR                default top distribution directory (build artifacts)
#     CND_BUILDDIR               default top build directory (object files, ...)
#     CONF                       name of current configuration
#     CND_PLATFORM_${CONF}       platform name (current configuration)
#     CND_ARTIFACT_DIR_${CONF}   directory of build artifact (current configuration)
#     CND_ARTIFACT_NAME_${CONF}  name of build artifact (current configuration)
#     CND_ARTIFACT_PATH_${CONF}  path to build artifact (current configuration)
#     CND_PACKAGE_DIR_${CONF}    directory of package (current configuration)
#     CND_PACKAGE_NAME_${CONF}   name of package (current configuration)
#     CND_PACKAGE_PATH_${CONF}   path to package (current configuration)
#
# NOCDDL


# Environment 
MKDIR=mkdir
CP=cp
CCADMIN=CCadmin

# Misc directories
LIB_DIR=$(CND_BASEDIR)/dist/lib/
SRC_DIR=$(CND_BASEDIR)/src/
PYTHON_PACKAGE_DIR=$(SRC_DIR)/python/molpher/swig_wrappers
TBB_LIBS=$(CND_BASEDIR)/deps/tbb/lib/intel64/gcc4.4/

# filenames
LIBNAME="libmolpher.so"

# run swig
RUN_SWIG=$(or $(SWIG),$(findstring SWIG,$(CONF)))
IS_RELEASE=$(findstring Release,$(CONF))

# build
build: .build-post

.build-pre:
# Add your pre 'build' code here...
	mkdir -p $(LIB_DIR)
ifneq (,$(RUN_SWIG))
	swig3.0 -python -Iinclude/ -Wall -c++ -outdir $(PYTHON_PACKAGE_DIR) -o $(SRC_DIR)/swig/molpher_wrap.cpp $(CND_BASEDIR)/include/molpher.i
endif
	
.build-post: .build-impl
# Add your post 'build' code here...

# copy library files to the distribution directory for testing purposes
ifneq (,$(IS_RELEASE))	
	cp $(TBB_LIBS)libtbb.so.2 $(LIB_DIR)libtbb.so.2
	ln -sf $(LIB_DIR)libtbb.so.2 $(LIB_DIR)libtbb.so
	cp $(TBB_LIBS)libtbbmalloc.so.2 $(LIB_DIR)libtbbmalloc.so.2
	ln -sf $(LIB_DIR)libtbbmalloc.so.2 $(LIB_DIR)libtbbmalloc.so
	cp $(TBB_LIBS)libtbbmalloc_proxy.so.2 $(LIB_DIR)libtbbmalloc_proxy.so.2
	ln -sf $(LIB_DIR)libtbbmalloc_proxy.so.2 $(LIB_DIR)libtbbmalloc_proxy.so
else
	cp $(TBB_LIBS)libtbb_debug.so.2 $(LIB_DIR)libtbb_debug.so.2
	ln -sf $(LIB_DIR)libtbb_debug.so.2 $(LIB_DIR)libtbb_debug.so
	cp $(TBB_LIBS)libtbbmalloc_debug.so.2 $(LIB_DIR)libtbbmalloc_debug.so.2
	ln -sf $(LIB_DIR)libtbbmalloc_debug.so.2 $(LIB_DIR)libtbbmalloc_debug.so
	cp $(TBB_LIBS)libtbbmalloc_proxy_debug.so.2 $(LIB_DIR)libtbbmalloc_proxy_debug.so.2
	ln -sf $(LIB_DIR)libtbbmalloc_proxy_debug.so.2 $(LIB_DIR)libtbbmalloc_proxy_debug.so
endif
	
	cp $(CND_BASEDIR)/res/SAScore.dat $(PYTHON_PACKAGE_DIR)/SAScore.dat

# clean
clean: .clean-post
	rm -rf $(LIB_DIR)
	rm -f $(PYTHON_PACKAGE_DIR)/*.dat
	rm -rf dist/
	rm -rf build/
	rm -f $(PYTHON_PACKAGE_DIR)/*.so

.clean-pre:
# Add your pre 'clean' code here...

.clean-post: .clean-impl
# Add your post 'clean' code here...


# clobber
clobber: .clobber-post

.clobber-pre:
# Add your pre 'clobber' code here...

.clobber-post: .clobber-impl
# Add your post 'clobber' code here...


# all
all: .all-post

.all-pre:
# Add your pre 'all' code here...

.all-post: .all-impl
# Add your post 'all' code here...


# build tests
build-tests: .build-tests-post

.build-tests-pre:
# Add your pre 'build-tests' code here...

.build-tests-post: .build-tests-impl
# Add your post 'build-tests' code here...


# run tests
test: .test-post

.test-pre: build-tests
# Add your pre 'test' code here...

.test-post: .test-impl
# Add your post 'test' code here...


# help
help: .help-post

.help-pre:
# Add your pre 'help' code here...

.help-post: .help-impl
# Add your post 'help' code here...



# include project implementation makefile
include nbproject/Makefile-impl.mk

# include project make variables
include nbproject/Makefile-variables.mk
