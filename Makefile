all: myMinDS

myMinDS: myMinDS.cpp graph.h localSearch.h constants.h dsBuckets.h\
					myBijection.h myBuckets.h hugeInt.h dominationHash.h
	g++ -std=gnu++0x -O3 -static myMinDS.cpp -o myMinDS

clean: rm -f *~
