#include <iostream>
#include <stdlib.h>
#include <iomanip>

#include <spinlib/fptype.h>

using namespace std;

void angleToVector(const double phi, const double theta)
{
//   cout << phi << '\t' << theta << endl;
  cout << setw(10) << sin(theta) * cos(phi) << '\t' << setw(10) << sin(theta) * sin(phi) << '\t' << setw(10) << cos(theta) << endl;
}

int main(int argc, char* argv[]) {
  const int knots = argc > 1 ? atoi(argv[1]) : 5;

  double piHalf = M_PI / 2.0;

  for (int k = 0; k < knots; ++k) {
    fp theta = piHalf * fp(k) / fp(knots - 1);
    if (k == 0) {
      angleToVector(0, theta);
    } else {
      for (int q = 0; q < k; ++q) {
        fp phi = piHalf * fp(q) / fp(k);
        for(int o = 0; o < 4; ++o) {
          angleToVector(phi + fp(o) * piHalf, theta);
        }
      }
    }
  }

  return 0;
}
