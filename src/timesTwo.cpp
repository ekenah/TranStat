#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp:export]
int timeTwo(int x) {
	return x * 2;
}
