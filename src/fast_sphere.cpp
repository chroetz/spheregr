#include <Rcpp.h>
using namespace Rcpp;

//' Converts a 2D tangentspace vector representation to 3D
//'
//' @param x vetor length 3
//' @param a vetor length 2
//' @return vector length 3
// [[Rcpp::export]]
NumericVector get_v_rcpp(NumericVector x, NumericVector a) {
  double b = 0;
  if (abs(x[0]-1) + abs(x[1]) + abs(x[2]) < 1e-7) b = 1;

  double u2[3];
  double u2p[3];
  double c = x[0] + b*x[1];
  u2p[0] = 1 - c*x[0];
  u2p[1] = b - c*x[1];
  u2p[2] = -c*x[2];

  double nrm = sqrt(u2p[0]*u2p[0] + u2p[1]*u2p[1]+ u2p[2]*u2p[2]);
  for (int i = 0; i < 3; ++i) u2[i] = u2p[i] / nrm;

  double u3[3];
  double u3p[3];

  u3p[0] = -x[2]*x[0] - u2[2]*u2[0];
  u3p[1] = -x[2]*x[1] - u2[2]*u2[1];
  u3p[2] = 1 - x[2]*x[2] - u2[2]*u2[2];
  nrm = sqrt(u3p[0]*u3p[0] + u3p[1]*u3p[1]+ u3p[2]*u3p[2]);
  if (nrm < 1e-7) {
    u3p[0] = -x[1]*x[0] - u2[1]*u2[0];
    u3p[1] = 1 -x[1]*x[1] - u2[1]*u2[1];
    u3p[2] = x[1]*x[2] - u2[1]*u2[2];
    nrm = sqrt(u3p[0]*u3p[0] + u3p[1]*u3p[1]+ u3p[2]*u3p[2]);
  }
  for (int i = 0; i < 3; ++i) u3[i] = u3p[i] / nrm;

  NumericVector out(3);
  for (int i=0; i<3; ++i) {
    out[i] = a[0]*u2[i] + a[1]*u3[i];
  }
  return out;
}



//' Convert angle representation to R3.
//'
//' @param x vetor length 2
//' @return vector length 3
// [[Rcpp::export]]
NumericVector angle2R3_1(NumericVector x) {
  NumericVector out(3);
  out[0] = sin(x[0])*cos(x[1]);
  out[1] = sin(x[0])*sin(x[1]);
  out[2] = cos(x[0]);
  return out;
}

