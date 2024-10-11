rem %1 - $(SolutionDir)
rem %2 - $(Configuration)

cd %1\libs

@echo "zlib-ng"
cd zlib-ng 
cmake -B build-vs/zlib-ng -S . -DZLIB_COMPAT=ON
cmake --build build-vs/zlib-ng --config %2
cd ..

@echo "isa-l"
cd isa-l
nmake -f Makefile.nmake
cd ..

@echo "libdeflate"
cd libdeflate
cmake -B build
cmake --build build --config %2
cd ..
