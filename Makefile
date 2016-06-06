DLLS := fftw32.dll fftw64.dll lapack32.dll lapack64.dll

math.zip: math.dyalog $(DLLS)
	zip $@ $^

fftw32.dll: | build-fftw32
	cd build-fftw32 && ../fftw/configure --host=mingw32 --disable-alloca --with-our-malloc16 --with-windows-f77-mangling --enable-shared --disable-static --with-incoming-stack-boundary=2 --enable-sse2 --enable-avx CC=i686-w64-mingw32-gcc && $(MAKE) CCLD="i686-w64-mingw32-gcc -Wc,-static-libgcc"
	cp build-fftw32/.libs/libfftw3-3.dll $@

fftw64.dll: | build-fftw64
	cd build-fftw64 && ../fftw/configure --host=mingw32 --disable-alloca --with-our-malloc16 --with-windows-f77-mangling --enable-shared --disable-static --with-incoming-stack-boundary=2 --enable-sse2 --enable-avx CC=x86_64-w64-mingw32-gcc && $(MAKE) CCLD="x86_64-w64-mingw32-gcc -Wc,-static-libgcc"
	cp build-fftw64/.libs/libfftw3-3.dll $@

lapack32.dll: | build-lapack32
	$(MAKE) -C build-lapack32 -f ../Makefile.lapack CC=i686-w64-mingw32-gcc FC=i686-w64-mingw32-gfortran
	cp build-lapack32/lapack.dll $@

lapack64.dll: | build-lapack64
	$(MAKE) -C build-lapack64 -f ../Makefile.lapack CC=x86_64-w64-mingw32-gcc FC=x86_64-w64-mingw32-gfortran
	cp build-lapack64/lapack.dll $@

build-%:
	mkdir $@
