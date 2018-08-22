libdistributions.so: distributions_lib.cxx
	g++ `root-config --libs --cflags` -lMinuit -O3 -fPIC -c distributions_lib.cxx -o distributions.o
	g++ `root-config --libs --cflags` -lMinuit -O3 -shared -Wl,-soname,libdistributions.so -o libdistributions.so distributions.o

clean:
	rm *.o *.so

