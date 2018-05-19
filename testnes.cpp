#include <iostream>

#include "nanoes.hpp"


#include <limits>

#include <stdio.h>
#include <iostream>
#include <time.h>


int main() {
	nanoes::runtime rt(10000);
	rt.set("print", std::function<int(void*, std::string)>([](auto ths,auto f) {
		std::cout << "Print[" << f << "]\n";
		return 0;
	}));
	rt.set("xstr", std::string("global X value"));
	rt.set("true", true);
	rt.set("false", false);
#if 1
	// test ext-call and string-const
	rt.eval(R"-(
		print("Hello world");
	)-");
	// test string+bool add
	rt.eval(R"-(
		print("Hello world:"+true);
	)-");
	// test string+num add
	rt.eval(R"-(
		print("Hello world"+1);
	)-");
	// test string + subblock arith and operator precedence
	rt.eval(R"-(
		print("Should be 2:"+(1+1));
		print("Should be 2:"+(3-1));
		print("Should be 7:" +(1+2*3));
		print("Should be 11:"+(1+2*3+4));
		print("Should be 27:"+(1+2*3+4*5));
		print("Should be 33:"+(1+2*3+4*5+6));
	)-");
	// comparasion ops
	rt.eval(R"-(
		print("Smaller? Y:"+(1 < 2));
		print("Smaller?: N"+(1 < 1));
		print("LEQ? Y:"+(1 <= 1));
		print("LEQ? N:"+(2 <= 1));
	)-");
	// test IF code
	rt.eval(R"-(
		if (1<2) {
			print("smaller should run\n");
		} else {
			print("smaller should NOT run\n");
		}
		if (2<1) {
			print("smaller should NOT run\n");
		} else {
			print("smaller should run\n");
		}
	)-");
#endif

#if 1
	// identifier scoping location.
	rt.eval(R"-(
		dummy("ARG VALUE");
		function dummy(xstr) {
			print("want arg not global ? -> "+xstr);
		}
	)-");
	// finally time the fibonacchi
	double start = clock();
	rt.eval(R"-(
		print("fib 12:"+fib(12));
		function fib(x) {
			if (x<2) {
				return x;
			} else {
				return fib(x-2)+fib(x-1);
			}
		}
	)-");
	double end = clock();
	printf("Time taken: %f\n",(end-start)/CLOCKS_PER_SEC);
#endif

#if 0
	// Early GC test(pre-interning string constants will produce garbage)
	rt.eval(R"-(
		while (true) {
			print("Loooooop...\n");
		}
	)-");
#endif

	return 0;
}
