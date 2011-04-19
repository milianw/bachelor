#include <iostream>

#include <spinlib/fptype.h>

using namespace std;

int main() {
  const int knots = 5;

  double piHalf = M_PI / 2.0;

  for (int k = 0; k < knots; ++k) {
    fp theta = piHalf * fp(k) / fp(knots - 1);
    if (k == 0) {
      cout << 0 << '\t' << theta << endl;
    } else {
      for (int q = 0; q < k; ++q) {
        fp phi = piHalf * fp(q) / fp(k);
        for(int o = 0; o < 4; ++o) {
          cout << (phi + fp(o) * piHalf) << '\t' << theta << endl;
        }
      }
    }
  }

  return 0;
}
