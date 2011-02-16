#include <iostream>
#include <bitset>
#include <cmath>

using namespace std;

const int nprotons = 10;
const int dimension = pow(2, nprotons + 1);
bitset<nprotons+1>* m_states;

bool contributesBit(int i, int j, int k)
{
  for(int l = 0; l < nprotons; ++l) {
    if (l == k) {
      continue;
    }
    if (m_states[i][l] != m_states[j][l]) {
      return false;
    }
  }
  return true;
}

bool contributes(int i, int j, int k)
{
  // states are equal if: all bits except for k-bit are equal
  // k-bit == 2^k = 0001000
  //                   ^k = 4
  int kBit = (1 << k);
  if (true) {
    // ignore electron bit
    // electronBit == 2^nprotons == 100000...
    int electronBit = (1 << nprotons);
    // essentially a fast variant of:
    // i % electronBit
    i &= electronBit - 1;
    j &= electronBit - 1;
  }
  return (i | kBit) == (j | kBit);
}

bool state(int i, int k)
{
  // k-bit == 2^k = 0001000
  //                   ^k = 4
  const int kPow = (1 << k);
  return i & kPow;
}

int main() {
  m_states = new bitset<nprotons+1>[dimension];
  //This includes all spin half species (protons+1 electron)
  for (int i = 0; i < dimension; ++i) {
    //for each bit: 0=+0.5, 1=-0.5
    m_states[i] = i;
  }

  for(int i = 0; i < dimension; ++i) {
    for(int j = 0; j < dimension; ++j) {
      for(int k = 0; k < nprotons; ++k) {
        bool bit = contributesBit(i, j, k);
        bool alt = contributes(i, j, k);
        if (m_states[i][k] != state(i, k)) {
          cout << "state unequal:" << i << ", " << k << ":" << m_states[i][k] << " VS " << state(i, k) << endl;
          return 2;
        }
        if (m_states[j][k] != state(j, k)) {
          cout << "state unequal:" << j << ", " << k << ":" << m_states[j][k] << " VS " << state(j, k) << endl;
          return 2;
        }
        const char * msg;
        bool failed = false;
        if (bit == alt) {
          msg = "";
        } else {
          msg = "\tFAIL";
          failed = true;
        }
        if (failed) {
          cout << "i = " << i << " = " << m_states[i] << "\t, j =" << j << " = " << m_states[j] << "\t, k = " << k << "\t:\t" << bit << '\t' << alt << msg << endl;
          return 1;
        }
      }
    }
  }

  return 0;
}