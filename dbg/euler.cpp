#include <iostream>

#include <Eigen/SVD>
#include <Eigen/Dense>
#include <map>

using namespace std;
using namespace Eigen;

Matrix3d RotZ(double angle)
{
  return (Matrix3d() << cos(angle), -sin(angle), 0,
                       sin(angle), cos(angle), 0,
                       0, 0, 1).finished();
}

Matrix3d RotY(double angle)
{
  return (Matrix3d() << cos(angle), 0, sin(angle),
                        0, 1, 0,
                        -sin(angle), 0, cos(angle)).finished();
}

Matrix3d RotX(double angle)
{
  return (Matrix3d() << 1, 0, 0,
                        0, cos(angle), -sin(angle),
                        0, sin(angle), cos(angle)).finished();
}

const double degree = M_PI / 180;
const double radians = 180 / M_PI;

Vector3d eulerAngles(const Matrix3d& rot)
{
  double beta = acos(rot(2, 2));
  double sinBeta = sqrt(1 - rot(2, 2) * rot(2, 2));
  double alpha = asin(rot(2, 1) / sinBeta);
  if (isnan(alpha)) {
    alpha = acos(rot(2, 0) / sinBeta);
  }
  double gamma = asin(rot(1, 2) / sinBeta);
  if (isnan(gamma)) {
    gamma = acos(rot(0, 2) / sinBeta);
  }
  return (Vector3d() << alpha, beta, gamma).finished();
}

Matrix3cd eigenVectors(const Matrix3d& m)
{
  return EigenSolver<Matrix3d>(m).eigenvectors();
}

double angle(const Vector3d& l, const Vector3d& r)
{
  return acos((l.dot(r) / (l.norm() * r.norm())));
}

Vector3d findAngles(const Matrix3cd& ref, const Matrix3cd& other)
{
  Matrix3d rot;
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < 3; ++j) {
      const Vector3cd& l = ref.col(i);
      const Vector3cd& r = other.col(j);
//       cout << l.norm() << ", " << r.norm() << ": " << (l.norm() * r.norm()) << endl;
      rot(i, j) = (l.dot(r) / (l.norm() * r.norm())).real();
    }
  }
  cout << "constructed rot via dot angles:" << endl;
  const double det = rot.determinant();
  if (det < 0) {
    rot.col(0) = -rot.col(0);
  }
  cout << rot << endl << "det = " << rot.determinant() << endl;
  const Vector3d euler = eulerAngles(rot);
  cout << "euler angles:" << euler.transpose() * radians << endl << endl;
  return euler;
}

Vector3d tryAngles(const Matrix3cd& ref, const Matrix3cd& other, int i, int j, int k)
{
  // 0, 2, 1
  cout << endl << "## sorted: " << i << ", " << j << ", " << k << endl;
  Matrix3cd sorted;
  sorted.col(0) = other.col(i);
  sorted.col(1) = other.col(j);
  sorted.col(2) = other.col(k);
  cout << "sorted eigenvectors:" << endl << sorted << endl << endl;
  return findAngles(ref, sorted);
}

int main() {
  const Matrix3d base = (Matrix3d() <<
//                           1, 1, 6,
//                           1, 1, 7,
//                           6, 7, -2
                          0.3856,      -0.4518,      -0.0070,
                          -0.4518,      -0.0268,       0.0090,
                          -0.0070,       0.0090,      -0.3431
                        ).finished();

/*
  const double alpha = 30;
  const double beta = 70;
  const double gamma = 80;
  const Matrix3d rot = RotZ(alpha * degree) * RotY(beta * degree) * RotZ(gamma * degree);
  cout << "rotation matrix: R = " << endl << rot << endl;
  cout << "euler angles:" << eulerAngles(rot).transpose() * radians << endl;
  cout << "expected are:" << alpha << ", " << beta << ", " << gamma << endl;
  const Matrix3d rotated = rot * base * rot.transpose();
  cout << "rotated base matrix: M' = " << endl << rotated << endl;
  cout << endl;
  const Matrix3cd RV = eigenVectors(rotated);
  cout << "base eigenvectors, rotated:" << endl << (rot * BV) << endl;
  cout << "rotated eigenvectors:" << endl << RV << endl;
*/
  const Matrix3d refFrame = (Matrix3d() << 
                            2.0055600,   -0.0000785,    0.0009609,
                            -0.0000448,    2.0060101,    0.0008908,
                            0.0009585,    0.0008655,    2.0028950
                            ).finished();
  cout << "reference matrix:" << endl << refFrame << endl;
  const Matrix3cd RV = eigenVectors(refFrame);
  cout << "reference matrix eigenvectors:" << endl << RV << endl;
  const Vector3d evals = EigenSolver<Matrix3d>(refFrame).eigenvalues().real();
  cout << "eigen values of refFrame:" << evals.transpose() << endl;
  Matrix3cd RV_sorted;
  map<double, int> sortedKeys;
  sortedKeys[evals(0)] = 0;
  sortedKeys[evals(1)] = 1;
  sortedKeys[evals(2)] = 2;
  std::map< double, int >::iterator it = sortedKeys.begin();
  RV_sorted.col(0) = RV.col(it->second);
  RV_sorted.col(1) = RV.col((++it)->second);
  RV_sorted.col(2) = RV.col((++it)->second);
  cout << "sorted eigenvectors of reference:" << endl << RV_sorted << endl << endl;

  cout << "base matrix: M = " << endl << base << endl;
  const Matrix3cd BV = eigenVectors(base);
  cout << "base eigenvectors:" << endl << BV << endl;
  const Vector3d BE = EigenSolver<Matrix3d>(base).eigenvalues().real();
  cout << "eigen values of base: " << BE.transpose() << endl;

  int minIndizes[3] = {0, 1, 2};
  Vector3d bestEuler;
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < 3; ++j) {
      if (i == j) continue;
      for(int k = 0; k < 3; ++k) {
        if (k == i || k == j) continue;
        Vector3d euler = tryAngles(RV_sorted, BV, i, j, k);
        if (i == 0 && j == 1 && k == 2) {
          bestEuler = euler;
        } else if (euler(0) < bestEuler(0)) {
          bestEuler = euler;
          minIndizes[0] = i;
          minIndizes[1] = j;
          minIndizes[2] = k;
        }
      }
    }
  }

  cout << "euler angles with orca algorithm: " << bestEuler.transpose() * radians << endl;
  cout << "eigen values: " << BE(minIndizes[0]) << ", " << BE(minIndizes[1]) << ", " << BE(minIndizes[2]) << endl;

  return 1;

  
  /*
  
  
  
  cout << "now lets try to find out the rotation matrix from base and rotated base:" << endl;
  const Matrix3d T = base.inverse() * rotated;
  cout << "transformation: T = M.transpose() * M' = " << endl << T << endl;
  cout << "applied transformation: M' = M * T = " << endl << base * T << endl;
  cout << endl;

  const Matrix3d Trot = Affine3d(T).rotation();
  cout << "rotation of T: U = " << endl << Trot << endl;
  cout << "applied rotation: M' = U * M * U.transpose() = " << endl << Trot * base * Trot.transpose() << endl;
  cout << endl;
  cout << "euler angles:" << endl << eulerAngles(Trot).transpose() * radians << endl;
  cout << "expected are:" << alpha << ", " << beta << ", " << gamma << endl;
  cout << endl << endl;


  cout << "perpendicular? x y " << angle(BV.col(0).real(), BV.col(1).real()) * radians << endl;
  cout << "perpendicular? y z " << angle(BV.col(1).real(), BV.col(2).real()) * radians << endl;
  cout << "perpendicular? x z " << angle(BV.col(0).real(), BV.col(2).real()) * radians << endl;
  cout << endl;

  cout << "z-axis angle beta =" << angle(base.col(0).real(), rotated.col(0).real()) * radians << endl;

  cout << angle((RotY(beta * degree) * base * RotY(beta * degree).transpose()).col(0).real(), rotated.col(0).real()) * radians << endl;

  {
  Matrix3d tmp = base;
  const double beta = eulerAngles(eigenVectors(tmp).real())(1) + eulerAngles(RV.real())(1);
  cout << "beta = " << beta << '\t' << beta * radians << endl;
  tmp = RotZ(beta) * tmp * RotZ(beta).transpose();
  cout << eulerAngles(eigenVectors(tmp).real()).transpose() * radians << endl;
  cout << tmp << endl;
  }
  */
  return 1;
}
