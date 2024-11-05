all: bkc bkc_dump

dummy := $(shell git submodule update --init --recursive)

BKC_MAIN_DIR = src/bkc
BKC_DUMP_DIR = src/bkc_dump
BKC_COMMON_DIR = src/common
BKC_LIBS_DIR = libs
LIBS_DIR = . #/usr/local/lib
INCLUDE_DIR= libs
ZLIB_INCLUDE_DIR = libs/zlib-ng/build-g++
ZLIB_INCLUDE_DIR_FOR_FILE_WRAPPER = libs/zlib-ng/build-g++/zlib-ng
ZSTD_INCLUDE_DIR = libs/zstd/lib
LIBDEFLATE_INCLUDE_DIR = libs/libdeflate
ISAL_INCLUDE_DIR = libs/isa-l/include
MIMALLOC_INLUCDE_DIR = libs/mimalloc/include

SHARED_INCLUDE_DIR ?= ./shared/

BKC_OUT_BIN_DIR ?= bin

MIMALLOC_OBJ=libs/mimalloc/mimalloc.o

ifdef MSVC     # Avoid the MingW/Cygwin sections
    UNAME_S := Windows
else                          # If uname not available => 'not'
    UNAME_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
    UNAME_M := $(shell uname -m)
endif

D_OS =
D_ARCH =

ifeq ($(UNAME_S),Darwin)
	D_OS=MACOS
	ifeq ($(UNAME_M),arm64)
		D_ARCH=ARM64
	else
		D_ARCH=X64
	endif
else
	D_OS=LINUX
	D_ARCH=X64

	ifeq ($(UNAME_M),arm64)
		D_ARCH=ARM64
	endif
	ifeq ($(UNAME_M),aarch64)
		D_ARCH=ARM64
	endif
endif

#ifeq ($(D_ARCH),X64)
#	dummy_install_nasm := $(shell \
#	if [ ! -f build_tools/nasm/nasm ]; then \
#		cd build_tools/nasm && ./autogen.sh && ./configure && make -j; \
#	fi)
#endif

#NASM_V := $(shell build_tools/nasm/nasm --version 2>/dev/null)

CMAKE_OSX_SYSROOT_FLAG =
ifeq ($(UNAME_S),Darwin)
	SDK_PATH := $(shell $(CXX) -v 2>&1 | grep -- '--with-sysroot' | sed -E 's/.*--with-sysroot=([^ ]+).*/\1/')
	CMAKE_OSX_SYSROOT_FLAG := -DCMAKE_OSX_SYSROOT=$(SDK_PATH)
endif

CPU_FLAGS =
STATIC_LFLAGS =
PLATFORM_SPECIFIC_FLAGS =

#in some cases we can have different results on ARM
#I guess this is exactly the same as here: https://bugs.mysql.com/bug.php?id=82760
ifeq ($(D_ARCH),ARM64)
	PLATFORM_SPECIFIC_FLAGS = -ffp-contract=off
endif

ifeq ($(D_OS),MACOS)
	ifeq ($(D_ARCH),ARM64)
		CPU_FLAGS = -march=armv8.4-a
	else
		CPU_FLAGS = -m64
	endif
	STATIC_LFLAGS = -static-libgcc -static-libstdc++ -pthread
else
	ifeq ($(D_ARCH),ARM64)
		CPU_FLAGS = -march=armv8-a
		STATIC_LFLAGS = -static-libgcc -static-libstdc++ -lpthread
	else
		CPU_FLAGS = -m64
		STATIC_LFLAGS = -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
	endif
endif

LIB_ZLIB=$(BKC_LIBS_DIR)/zlib-ng/build-g++/zlib-ng/libz.a
LIB_ISAL=$(BKC_LIBS_DIR)/isa-l/bin/isa-l.a
LIB_ZSTD=$(BKC_LIBS_DIR)/zstd/lib/libzstd.a
LIB_LIBDEFLATE=$(BKC_LIBS_DIR)/libdeflate/build/libdeflate.a

LIB_GZ=$(LIB_ZLIB)

REFRESH_FLAGS = 
ifeq ($(UNAME_S),Linux)
	ifeq ($(UNAME_M),x86_64)
		REFRESH_FLAGS +=-DARCH_X64
#		ifdef NASM_V
#			LIB_GZ=$(LIB_ISAL)
#			REFRESH_FLAGS += -DREFRESH_USE_IGZIP -I $(ISAL_INCLUDE_DIR)
#		else
			REFRESH_FLAGS +=-DREFRESH_USE_ZLIB
#		endif
	else
		REFRESH_FLAGS +=-DREFRESH_USE_ZLIB
	endif
else
	REFRESH_FLAGS +=-DREFRESH_USE_ZLIB
endif

CFLAGS	= -fPIC -Wall -O3 $(PLATFORM_SPECIFIC_FLAGS) $(CPU_FLAGS) $(REFRESH_FLAGS) -std=c++20 -pthread -I $(SHARED_INCLUDE_DIR) -I $(ZLIB_INCLUDE_DIR) -I $(ZLIB_INCLUDE_DIR_FOR_FILE_WRAPPER) -I $(INCLUDE_DIR) -I $(MIMALLOC_INLUCDE_DIR) -I $(ZSTD_INCLUDE_DIR) -I $(LIBDEFLATE_INCLUDE_DIR) -fpermissive
CLINK	= -lm -std=c++20 -lpthread

release: CLINK = -lm -std=c++20 $(STATIC_LFLAGS)

release: CFLAGS	= -fPIC -Wall -O3 -DNDEBUG $(PLATFORM_SPECIFIC_FLAGS) $(CPU_FLAGS) $(REFRESH_FLAGS) -std=c++20 -pthread -I $(SHARED_INCLUDE_DIR) -I $(ZLIB_INCLUDE_DIR) -I $(ZLIB_INCLUDE_DIR_FOR_FILE_WRAPPER) -I $(INCLUDE_DIR) -I $(MIMALLOC_INLUCDE_DIR) -I $(ZSTD_INCLUDE_DIR) -fpermissive
#release: satc_undump satc_filter satc_to_fasta fafq_filter
release: all

debug: CFLAGS	= -fPIC -Wall -O0 -g $(PLATFORM_SPECIFIC_FLAGS) $(CPU_FLAGS) $(REFRESH_FLAGS) -std=c++20 -pthread -I $(SHARED_INCLUDE_DIR) -I $(ZLIB_INCLUDE_DIR) -I $(ZLIB_INCLUDE_DIR_FOR_FILE_WRAPPER) -I $(INCLUDE_DIR) -I $(MIMALLOC_INLUCDE_DIR) -I $(ZSTD_INCLUDE_DIR) -fpermissive
debug: all

ifeq ($(UNAME_S),Linux)
	CLINK+=-fabi-version=6
endif


# default install location (binary placed in the /bin folder)
prefix      = /usr/local

# optional install location
exec_prefix = $(prefix)

$(LIB_ZLIB):
	cd $(BKC_LIBS_DIR)/zlib-ng; cmake $(CMAKE_OSX_SYSROOT_FLAG) -DCMAKE_CXX_COMPILER=$(CXX) -DCMAKE_C_COMPILER=$(CC) -B build-g++/zlib-ng -S . -DZLIB_COMPAT=ON; cmake --build build-g++/zlib-ng --config Release --target zlibstatic

#$(LIB_ISAL):
#	cd $(BKC_LIBS_DIR)/isa-l && PATH=../../build_tools/nasm:$$PATH make -f Makefile.unx

$(LIB_ZSTD):
	cd $(BKC_LIBS_DIR)/zstd; make -j
	
$(LIB_LIBDEFLATE):
	cd $(BKC_LIBS_DIR)/libdeflate; cmake $(CMAKE_OSX_SYSROOT_FLAG) -DCMAKE_CXX_COMPILER=$(CXX) -DCMAKE_C_COMPILER=$(CC) -B build && cmake --build build

$(MIMALLOC_OBJ):
	$(CXX) -DMI_MALLOC_OVERRIDE -O3 -DNDEBUG -fPIC -Wall -Wextra -Wno-unknown-pragmas -fvisibility=hidden -ftls-model=initial-exec -fno-builtin-malloc -c -I libs/mimalloc/include libs/mimalloc/src/static.c -o $(MIMALLOC_OBJ)

%.o: %.cpp $(LIB_ZLIB) $(LIB_GZ)
	$(CXX) $(CFLAGS) -c $< -o $@

bkc: $(BKC_OUT_BIN_DIR)/bkc

#seems to be important to have mimalloc first in link list, else segfault
$(BKC_OUT_BIN_DIR)/bkc: $(BKC_MAIN_DIR)/bkc.o \
	$(BKC_COMMON_DIR)/bkc_file.o \
	$(BKC_MAIN_DIR)/kmer_counter.o \
	$(BKC_MAIN_DIR)/fq_reader.o \
	$(BKC_MAIN_DIR)/memory_pool.o \
	$(BKC_COMMON_DIR)/utils.o \
	$(LIB_ZLIB) \
	$(LIB_ZSTD) \
	$(MIMALLOC_OBJ)
	-mkdir -p $(BKC_OUT_BIN_DIR)
	$(CXX) -o $@ \
	$(MIMALLOC_OBJ) \
	$(BKC_MAIN_DIR)/bkc.o \
	$(BKC_COMMON_DIR)/bkc_file.o \
	$(BKC_MAIN_DIR)/kmer_counter.o \
	$(BKC_MAIN_DIR)/fq_reader.o \
	$(BKC_MAIN_DIR)/memory_pool.o \
	$(BKC_COMMON_DIR)/utils.o \
	$(LIB_ZLIB) \
	$(LIB_ZSTD) \
	$(CLINK)

bkc_dump: $(BKC_OUT_BIN_DIR)/bkc_dump

#seems to be important to have mimalloc first in link list, else segfault
$(BKC_OUT_BIN_DIR)/bkc_dump: $(BKC_DUMP_DIR)/bkc_dump.o \
	$(BKC_COMMON_DIR)/bkc_file.o \
	$(BKC_DUMP_DIR)/dumper.o \
	$(LIB_ZSTD) \
	$(MIMALLOC_OBJ)
	-mkdir -p $(BKC_OUT_BIN_DIR)
	$(CXX) -o $@ \
	$(MIMALLOC_OBJ) \
	$(BKC_DUMP_DIR)/bkc_dump.o \
	$(BKC_COMMON_DIR)/bkc_file.o \
	$(BKC_DUMP_DIR)/dumper.o \
	$(LIB_ZSTD) \
	$(CLINK)

clean:
	-rm $(BKC_MAIN_DIR)/*.o
	-rm $(BKC_DUMP_DIR)/*.o
	-rm $(BKC_COMMON_DIR)/*.o
	-rm -rf $(BKC_OUT_BIN_DIR)
	cd $(BKC_LIBS_DIR)/zlib-ng && $(MAKE) -f Makefile.in clean
	cd $(BKC_LIBS_DIR)/zstd && make clean

