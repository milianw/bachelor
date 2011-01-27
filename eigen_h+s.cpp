#include <iostream>
#include <string>
#include <cmath>
#include <bitset>
#include <gsl/gsl_const_mksa.h>

#include <complex>
#include <Eigen/Dense>

#include <QtCore/QCoreApplication>
#include <QtCore/QStringList>

using namespace Eigen;
using namespace std;

//Physical constants, all in MKS units
const double hbar = GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR;
const double h = GSL_CONST_MKSA_PLANCKS_CONSTANT_H;
const double NUC_MAGNETON = GSL_CONST_MKSA_NUCLEAR_MAGNETON; // (J/T)
const double g_1H = 5.585694701 ;            // proton g factor
const double g_14N = 0.403761;               // N14 g factor
const double GAMMA_1H = 267.513E6 ;          // (rad/s.T)
const double GAMMA_14N =  19.331E6 ;         // (rad/s.T)
const double Bohrm = 9.27400949E-24;         // Bohr magneton in J/T
// Reminder: GAMMA_1H * hbar =  NUC_MAGNETON * g_1H
//const double B = 1;                          // Static magnetic field (T)
//const double B2 = ;                        // RF field
//const double LARMOR_1H = B * GAMMA_1H /2/M_PI / 1.0E6;

// Pauli Matrices
typedef complex<double> c_double;
const Matrix2cd pauliX = (Matrix2cd() << 0, 0.5,
                                         0.5, 0).finished();
const Matrix2cd pauliY = (Matrix2cd() << 0, c_double(0, -0.5),
                                         c_double(0, 0.5), 0).finished();
const Matrix2cd pauliZ = (Matrix2cd() << 0.5, 0,
                                         0, -0.5).finished();

#ifndef PROTONS
const int nprotons = 7;
#else
const int nprotons = PROTONS;
#endif

const int dimension = pow(2, nprotons + 1);

namespace input {

//gtensor=============================
static inline Matrix3cd gTensor()
{
  return Matrix3cd::Identity() * 2.0022838;
}

//static field========================
static inline Vector3cd staticBField(const double B)
{
  const double B0_norm = B; //0.2839; //field in Tesla
  double B0[3] = {0, 0, 1}; //field direction
  const double B_temp = sqrt( B0[0]*B0[0] + B0[1]*B0[1] + B0[2]*B0[2] );

  for (int i = 0; i<3; i++) {
    B0[i] = B0[i]*B0_norm/B_temp;
  }

  //cout << endl << "B0:" << " " << B0[0] << '\t' << B0[1] << '\t' << B0[2] << endl;
  //cout.precision(5);
  //cout << "\nStatic Field (T): " << sqrt(B0[0]*B0[0]+B0[1]*B0[1]+B0[2]*B0[2]) << endl;

  return Vector3cd(B0[0], B0[1], B0[2]);
}

//arbitrary A-tensor==================
static inline Matrix3cd aTensor()
{
  return Matrix3cd::Identity() * 1420;  //proton hyperfine coupling in T
}

} // namespace input

//BEGIN SpinHamiltonian

/**
 * TODO: check whether inlining some parts is noticeable
 *
 * notes on porting:
 *
 * ~~~~~~~~~~~
 * gsl_blas_zdotu (x, y, dotu)
 * => eigen: dotu = x.conjugate().dot(y)
 *
 * from eigen docs about .dot():
 * Note: If the scalar type is complex numbers, then this function returns the hermitian (sesquilinear) dot product,
 *       conjugate-linear in the first variable and linear in the second variable.
 *       but gsl_blas_zdotu seems to differ from this, hence use a.conjugate().dot(b) instead of a.dot(b)
 * NOTE: <orzel> in this very specific case, it might be that (a.transpose()*b)(0,0) is faster (taking the only element of the 1x1 matrix x^T.y
 * ~~~~~~~~~~~
 * gsl_blas_zgemv(CblasTrans, gsl_complex_rect(1,0), atensor, I, gsl_complex_rect(0,0), atensor_I);
 * => eigen: atensor_I = atensor.transpose() * I
 * ~~~~~~~~~~~
 * gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, gsl_complex_rect(1,0),
 *                eigenvectors, moments, gsl_complex_rect(1,0), intermediate);
 * docs: C = \alpha op(A) op(B) + \beta C
 * => eigen: intermediate = eigenvectors.adjoint() * moments + intermediate
 */
class SpinHamiltonian {
  public:
    /// @p B magnetic field in Tesla
    SpinHamiltonian(const double B);

    void calculate() const;

  private:
    /// nuclear Zeeman Hamiltonian component
    MatrixXcd nuclearZeeman() const;
    /// hyper fine Hamiltonian component
    MatrixXcd hyperFine() const;
    /// electron Zeeman Hamiltonian component
    MatrixXcd electronZeeman() const;

    /// return spin vector from pauli matrices
    inline Vector3cd spinVector(int i, int j, int k) const;

    /// all bits for states @p i, @p j, must match for nucleus @p k
    /// otherwise the integral will be zero anyways and we can skip
    /// <  i  |H|  j  >
    /// <10101|H|10001> = <1|H|0><1|1><0|0><0|0><1|1> => must be calculated
    ///    ^=k     ^=k
    /// <10100|H|10001> = <1|H|0><1|1><0|0><0|0><0|1>
    ///    ^=k     ^=k                            ^=0 => can be skipped
    inline bool stateContributes(const int i, const int j, const int k) const;

    /// moments
    /// TODO: better document
    MatrixXcd magneticMoments() const;

    inline c_double magneticMoment(const int i, const int j) const;

    /// probability matrix with coefficients (i, j) = |< psi_j | M | psi_i>|^2
    /// psi_i being the i-th eigen vector
    /// M being the magnetic moment matrix
    MatrixXd probabilityMatrix(const MatrixXcd& eigenVectors) const;

    const double m_B;
    Matrix3cd m_gTensor;
    Matrix3cd m_aTensor;
    Vector3cd m_staticBField;
    bitset<(nprotons+1)>* m_states;
};

SpinHamiltonian::SpinHamiltonian(const double B)
: m_B(B)
, m_gTensor(input::gTensor())
, m_aTensor(input::aTensor())
, m_staticBField(input::staticBField(B))
, m_states(new bitset<(nprotons+1)>[dimension])
{
  //This includes all spin half species (protons+1 electron)
  for (int i = 0; i < dimension; ++i) {
    ///TODO: what does this comment mean?
    //for each bit: 0=+0.5, 1=-0.5
    m_states[i] = i;
  }
}

inline Vector3cd SpinHamiltonian::spinVector(int i, int j, int k) const
{
  const int a = m_states[i][k]; //spin state of state k in row i
  const int b = m_states[j][k];  //spin state of state k in column j
  return (Vector3cd() << pauliX(a, b), pauliY(a, b), pauliZ(a, b)).finished();
}

inline bool SpinHamiltonian::stateContributes(const int i, const int j, const int k) const
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


MatrixXcd SpinHamiltonian::nuclearZeeman() const
{
  //Compute nZeeman============================================================  
  MatrixXcd nZeeman(dimension, dimension);
  nZeeman.setZero();
  //to turn off: return nZeeman;

  for (int i = 0; i < dimension; ++i) {
    for (int j = 0; j < dimension; ++j) {
      //nprotons is always the index of the electronic spin state

      if (m_states[i][nprotons] != m_states[j][nprotons]) {
        ///TODO: understand
        continue;  //matrix elements between different e states are zero
      }

      for (int k = 0; k < nprotons; ++k) {
        if (!stateContributes(i, j, k)) {
          continue;
        }

        // set cell to dot product of I and B
        ///TODO: .dot() or .conjugate().dot()
//         nZeeman(i, j) += m_staticBField.dot(spinVector(i, j, k));
        nZeeman(i, j) += m_staticBField.conjugate().dot(spinVector(i, j, k));
      }
    }
  }

  nZeeman *= (-1.0*g_1H*NUC_MAGNETON);

  // DEBUG:
  // cout << nZeeman << endl;

  return nZeeman;
}

MatrixXcd SpinHamiltonian::hyperFine() const
{
  //Compute Hyperfine couplings matrix=========================================
  MatrixXcd hyperfine(dimension, dimension);
  hyperfine.setZero();
  ///NOTE: why transpose?
  const MatrixXcd aTensorTransposed = m_aTensor.transpose();

  for (int i = 0; i < dimension; ++i) {
    for (int j = 0; j < dimension; ++j) {
      //compute elements of s vector
      const Vector3cd s = spinVector(i, j, nprotons);

      for (int k = 0; k < nprotons; ++k) {    //loop over nuclei
        if (!stateContributes(i, j, k)) {
          continue;
        }

        //skip this for now, all nuclei assigned the free proton isotropic coupling
        //multiply atensor by I
        //multiply s by atensor_I
        ///TODO: a.dot(b) or a.conjugate().dot(s)
        hyperfine(i, j) += (aTensorTransposed * spinVector(i, j, k)).conjugate().dot(s);
      }
    }
  }
  hyperfine *= h * 1.0E6;

  // DEBUG:
  // cout << hyperfine << endl;

  return hyperfine;
}

MatrixXcd SpinHamiltonian::electronZeeman() const
{
  //Compute eZeeman============================================================  
  MatrixXcd eZeeman(dimension, dimension);
  eZeeman.setZero();
  //first multiply the g tensor by the field
  const Vector3cd Bstaticg = m_gTensor.transpose() * m_staticBField;

  //depending on the convention, i might have to tranpose the gtensor here 
  //cout << "Bstaticg: \n";
  //for(int i=0; i<3; i++)
  //cout << GSL_REAL(gsl_vector_complex_get(Bstaticg, i)) << '\t';
  //cout << endl;
  for (int i = 0; i < dimension; ++i) {
    //cout << m_states[i] << " |";
    for (int j = 0; j < dimension; ++j) {
      //nprotons is always the index of the electron spin
      const int a = m_states[i][nprotons];
      const int b = m_states[j][nprotons];
      //cout << '\t' << a << b;
      if (!stateContributes(i, j, nprotons)) {
        continue;
      }

      eZeeman(i,j) += Bstaticg.transpose().dot(spinVector(i, j, nprotons));
    }
    //cout << endl;
  }
  eZeeman *= Bohrm;

  // DEBUG
  // cout << eZeeman << endl;

  return eZeeman;
}

MatrixXcd SpinHamiltonian::magneticMoments() const
{
  MatrixXcd moments(dimension, dimension);
  moments.setZero();

  for (int i = 0; i < dimension; ++i) {
    for (int j = 0; j < dimension; ++j) {
      //nprotons is always the index of the electronic spin state
      moments(i, j) = magneticMoment(i, j);

//       cout << i << '\t' << j << '\t' << "FINAL:" << '\t' << moments(i, j) << endl;
    }
  }
  return moments;
}

inline c_double SpinHamiltonian::magneticMoment(const int i, const int j) const
{
  c_double ret = 0;
  for (int k = 0; k < nprotons+1; ++k) {
    bool contributes = true;
    for (int l = 0; l < nprotons+1; ++l) {
      if (l==k) {
        continue;
      }
      if (m_states[i][l] != m_states[j][l]) {
        contributes = false;
        break;
      }
    }
    if (!contributes) {
      continue;
    }

    const int a = m_states[i][k];  //spin state of state k in row i
    const int b = m_states[j][k];  //spin state of state k in column j
    c_double xMoment = pauliX(a, b);

    if (k != nprotons) {
      xMoment *= -1.0 * g_1H * NUC_MAGNETON;
    } else {
      xMoment *= 2.023 * Bohrm;
    }

//         cout << i << '\t' << j << '\t' << k << '\t' << xMoment << endl;
    ret += xMoment;
  }
  return ret;
}

MatrixXd SpinHamiltonian::probabilityMatrix(const MatrixXcd& eigenVectors) const {
  const MatrixXcd moments = magneticMoments();
  //cout << moments << endl;
  MatrixXd probabilities(dimension, dimension);
  for(int i = 0; i < dimension; ++i) {
    for(int j = 0; j < dimension; ++j) {
      if (i != j) {
        probabilities(i, j) = ((eigenVectors.col(j).adjoint() * (moments * eigenVectors.col(i)))).norm();
      } else {
        probabilities(i, j) = 0;
      }
    }
  }
  return probabilities /= probabilities.maxCoeff();
}

void SpinHamiltonian::calculate() const
{
  //cout << "g = " << m_aTensor << endl;
  //cout << "B = " << m_staticBField << endl;
  //cout << "A = " << m_aTensor << endl;
  //cout << "Pauli Matrices:" << endl << pauliX << endl << pauliY << endl << pauliZ << endl;
  //cout << "Spin states: \n";
  //for(int i = 0; i < dimension; ++i) {
  //  cout << i+1 << '\t' << m_states[i] << endl;
  //}

  const MatrixXcd Hamiltonian = nuclearZeeman() + hyperFine() + electronZeeman();
  //cout << endl << endl << "~~~~~~~~~~~~~ Hamiltonian:" << endl;
  //cout << Hamiltonian << endl << "~~~~~~~~~~~~~" << endl;

  //Diagonalize the total Hamiltonian matrix===================================
  SelfAdjointEigenSolver<MatrixXcd> eigenSolver(Hamiltonian);
  const VectorXd eigenValues = eigenSolver.eigenvalues();
  const MatrixXcd eigenVectors = eigenSolver.eigenvectors();
  //cout << "eigenvectors:\n" << eigenVectors << endl;
  //cout << "\neigenvalues:\n" << eigenValues << endl;

  const MatrixXd probabilities = probabilityMatrix(eigenVectors);
  //cout << probabilities << endl;

  cout << "\n\nAllowed Transitions: \t\t\tB= " << fixed << m_B << endl;
  cout << "-------------------- " << endl;
  cout << "Transition\tFrequency (GHz)\t\tProbability" << endl;
  cout << "---------------------------------------------------\n";
  int transitions = 0;
  for (int i = 0;i < dimension; ++i) {
    for (int j = i + 1; j < dimension; ++j) {
      const double probability = probabilities(i, j);
      if (probability > 1.0E-6) {
          cout << fixed << i+1 << " -> " << j+1 << "\t\t";
          cout.precision(5);cout.width(10);
          // transition frequency:
          cout << right << (1.0/h/1.0E9 * abs(eigenValues(i) - eigenValues(j))) << "\t\t";
          cout.precision(8);
          cout << probability << endl;

          ++transitions;
        }
    }
  }
  cout << transitions << " transitions in total" << endl;
}

//END SpinHamiltonian

int main(int argc, char* argv[])
{
  //cout << 2.023 * Bohrm << endl;
  //cout << NUC_MAGNETON << endl;
  //cout << g_1H * NUC_MAGNETON << endl;
  //cout << GAMMA_1H * hbar << endl;

  /////////////////////////////////////////////////////////////////////////
  // The input file should be structured as follows:                 //
  // natoms n                     (number of coupled nuclei)         //
  // field x y z                  (The static field direction)       //
  // g                            (g tensor)                         //
  //   xx xy xz                              //
  //   yx yy yz                              //
  //   zx zy zz                              //
  // N                            (A tensor)                 //
  //   xx xy xz                              //
  //   yx yy yz                              //
  //   zx zy zz                              //
  // H                            (A tensor)                 //
  //   .........                             //
  //   .........                             //
  //   .........                             //
  // And so on...                                //
  /////////////////////////////////////////////////////////////////////////
  

  //Parsing input file
  //ifstream inputfile;
  //inputfile.open(argv[1]);
  // if(!inputfile)
  //   {
  //     cerr << "This program expects a single input file\n";
  //     cerr << "Check the comments in the source code for details\n"
  //     return(1);
  //   }


  //  for(int i=0;i<100;i++)
  //{

//     SpinHamiltonian h(0.362562);
//     const double B = 0.3;
    const double B = 0.3417757;
    SpinHamiltonian h(B);
    h.calculate();
  //}

  return 0;
}
