H5_INC_DIR = $(shell pwd)/hdf5/include
H5_LIB_DIR = $(shell pwd)/hdf5/lib

# ADDITIONAL_OPT_FLAGS = -march=znver4

CC = mpicc

OPT_FLAGS = -std=c11 -O3
DEBUG_FLAGS = -Wall -Wextra

LIB_FLAGS = -I$(H5_INC_DIR) -L$(H5_LIB_DIR) -lhdf5 -lz -lm

CFLAGS = $(OPT_FLAGS) $(ADDITIONAL_OPT_FLAGS) $(DEBUG_FLAGS) $(LIB_FLAGS)

SRC = $(wildcard src/*.c)
HEAD = $(wildcard include/*.h)
OBJ = $(patsubst src/%.c,obj/%.o,$(SRC))
HDF5_TAR = $(wildcard archives/hdf5*.tar.gz)

obj/%.o: src/%.c definitions.h $(HEAD)
	@mkdir -p obj
	@printf '\033[1;34mCC   %s\033[0m\n' "$<"
	@$(CC) -c $< -o $@ $(CFLAGS)

obj/main.o: main.c definitions.h $(HEAD)
	@mkdir -p obj
	@printf '\033[1;34mCC   %s\033[0m\n' "$<"
	@$(CC) -c $< -o $@ $(CFLAGS)

cfd: $(OBJ) obj/main.o
	@printf '\033[1;32mCCLD %s\033[0m\n' "$^"
	@$(CC) $^ -o $@  $(CFLAGS)

clean:
	rm -rf cfd obj

cleandata:
	rm -f *.h5

hdf5: $(HDF5_TAR)
	mkdir -p hdf5
	tar -xvzf $(HDF5_TAR) -C hdf5 --strip-components=2
	cd hdf5;\
	export ac_cv_lib_sz_SZ_BufftoBuffCompress=no;\
	export ac_cv_header_szlib_h=no;\
	CC=mpicc ./configure --enable-parallel --disable-shared --disable-szlib;\
	make;\
	make install;\
	cd ../;\
	mv hdf5/hdf5/include hdf5/;\
	mv hdf5/hdf5/lib hdf5/
