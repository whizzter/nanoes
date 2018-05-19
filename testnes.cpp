#include "nanoes.hpp"


#include <iostream>
#include <limits>

#include <stdio.h>
#include <iostream>
#include <time.h>


int main() {
	nanoes::runtime rt;
	rt.set("print", std::make_shared<std::function<int(void*, std::string)>>([](auto ths,auto f) {
		std::cout << "Print[" << f << "]\n";
		return 0;
	}));
	rt.set("xstr", std::string("global X value"));
	rt.set("true", true);
	rt.set("false", false);

	rt.eval(R"-(
		print("Hello world");
	)-");
	rt.eval(R"-(
		print("Hello world:"+true);
	)-");
	rt.eval(R"-(
		print("Hello world"+1);
	)-");
	rt.eval(R"-(
		print("Should be 2:"+(1+1));
		print("Should be 2:"+(3-1));
		print("Should be 7:" +(1+2*3));
		print("Should be 11:"+(1+2*3+4));
		print("Should be 27:"+(1+2*3+4*5));
		print("Should be 33:"+(1+2*3+4*5+6));
	)-");
	rt.eval(R"-(
		print("Smaller? Y:"+(1 < 2));
		print("Smaller?: N"+(1 < 1));
		print("LEQ? Y:"+(1 <= 1));
		print("LEQ? N:"+(2 <= 1));
	)-");
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
	return 0;
}
