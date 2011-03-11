/*
 * This file is part of my bachelor thesis.
 *
 * Copyright 2011 Milian Wolff <mail@milianw.de>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Library General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, write to the
 * Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef MW_BACHELOR_OPERATORSUM_P_H
#define MW_BACHELOR_OPERATORSUM_P_H

/**
 * Implementation of e.q. 2-77 on p 26/27
 *
 * TODO: extend for higher spins
 */

#define ASSERT(x) if(!(x)) { std::cerr << "assertion failed: " << #x << endl << " in " << __FUNCTION__ << " on line " << __LINE__ << endl; exit(1); }

class AbstractOperatorSummand
{
public:
  AbstractOperatorSummand() {}
  virtual ~AbstractOperatorSummand() {}

  virtual int dimension() const = 0;
  virtual c_fp operator()(int row, int col) const = 0;
  virtual AbstractOperatorSummand* copy() const = 0;
};

template<typename T>
class OperatorSummand : public AbstractOperatorSummand {
public:
  explicit OperatorSummand(const T& val)
  : m_val(val)
  { ASSERT(val.rows() == val.cols()); }
  virtual ~OperatorSummand()
  { }

  virtual int dimension() const
  {
    return m_val.rows();
  }

  virtual c_fp operator()(int row, int col) const
  {
    return m_val(row, col);
  }

  AbstractOperatorSummand* copy() const
  {
    return new OperatorSummand<T>(m_val);
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  const T m_val;
};

template<typename T>
const AbstractOperatorSummand* summandForType(const T& val)
{
  return new OperatorSummand<T>(val);
}

class OperatorSum : public AbstractOperatorSummand {
public:
  OperatorSum()
  : left(0), right(0)
  { }

  explicit OperatorSum(const AbstractOperatorSummand* l, const AbstractOperatorSummand* r)
  : left(l), right(r)
  { }

  template<typename LeftT, typename RightT>
  explicit OperatorSum(const LeftT& l, const RightT& r)
  : left(summandForType(l)), right(summandForType(r))
  { }

  OperatorSum(const OperatorSum& other)
  : left(other.left ? other.left->copy() : 0), right(other.right ? other.right->copy() : 0)
  { }

  virtual AbstractOperatorSummand* copy() const
  {
    return new OperatorSum(*this);
  }

  virtual ~OperatorSum()
  {
    delete left;
    delete right;
  }

  virtual int dimension() const
  {
    if (left && right) {
      return left->dimension() * right->dimension();
    } else if (left) {
      return left->dimension();
    } else if (right) {
      return right->dimension();
    } else {
      return 0;
    }
  }

  int rows() const
  {
    return dimension();
  }

  int cols() const
  {
    return dimension();
  }

  virtual c_fp operator()(int row, int col) const
  {
    ASSERT(row >= 0 && row < dimension());
    ASSERT(col >= 0 && col < dimension());

    if (!left) {
      return right->operator()(row, col);
    } else if (!right) {
      return left->operator()(row, col);
    }

    c_fp ret = 0;
    /// general notes
    // the kronecker product of A x B will result in a matrix with dim(A) submatrices
    // the submatrices are each defined by: A(i, j) * B
    // i, j can be obtained by row / dim(B) and col / dim(B) respectively
    const int rowA = row / right->dimension();
    const int colA = col / right->dimension();
    // analog for the index in B by replacing the division by a modulo operation
    const int rowB = row % left->dimension();
    const int colB = col % left->dimension();
    /// first term in operator sum:
    // J_1 x I = [[ J_1(0,0) * I , J_1(0, 1) * I, ..., J_1(0, n) * I ], ...]
    // J_1 is first operator, I is identity matrix of J_2's dimensionality
    // the kronecker product between J_1 and I result in a matrix
    // containing block diagonal sub matrices:
    // Each sub matrix is dim(I) x dim(I) sized and contains a single element
    // of J_1 on it's diagonal.
    // There are a total of dim(I) = dim(J_2) blocks in the matrix
    ///
    // check whether row, col lie on a diagonal of I, i.e. B
    if (rowB == colB) {
      // return value of A, with which the identity matrix gets scaled
      ret += left->operator()(rowA, colA);
    }
    /// second term in operator sum:
    // I x J_2 = [[ J, 0 ...], [0, J, ...], ... [..., J]]
    // J_2 is the second operator, I is the identity matrix of J_1's dimensionality
    // the kronecker product between J_2 and I result in a matrix
    // containing J_2 as submatrices in the diagonal (dim(I) times) and zero everywhere else
    ///
    // check whether row, col lie in on a diagonal of I, i.e. of A
    if (rowA == colA) {
      // return value of B (A scales by 1)
      ret += right->operator()(rowB, colB);
    }

    return ret;
  }

  template<typename T>
  OperatorSum operator+(const T& addee)
  {
    OperatorSum sum(*this);
    sum += addee;
    return sum;
  }

  template<typename T>
  OperatorSum& operator+=(const T& addee)
  {
    if (!left) {
      left = summandForType(addee);
    } else if (!right) {
      right = summandForType(addee);
    } else {
      right = new OperatorSum(right, summandForType(addee));
    }
    return *this;
  }

  MatrixXc toMatrix() const
  {
    MatrixXc ret(dimension(), dimension());
    for(int i = 0; i < dimension(); ++i) {
      for(int j = 0; j < dimension(); ++j) {
        ret(i, j) = this->operator()(i, j);
      }
    }
    return ret;
  }

  operator MatrixXc() const
  {
    return toMatrix();
  }

  const AbstractOperatorSummand* left;
  const AbstractOperatorSummand* right;

private:
  void operator=(const OperatorSum& other);
};

template<typename DataT>
class SpinOperator {
public:
  explicit SpinOperator()
  { }
  explicit SpinOperator(const DataT& _X, const DataT& _Y, const DataT& _Z)
  : X(_X), Y(_Y), Z(_Z)
  { }

  template<typename OtherDataT>
  SpinOperator& operator+=(const SpinOperator<OtherDataT>& other)
  {
    X += other.X;
    Y += other.Y;
    Z += other.Z;
  }

  /**
   * A * B, both being Spin Operators, equals
   * A.X (x) B.X + A.Y (x) B.Y + A.Z (x) B.Z
   *
   * with (x) denoting the direct product expansion
   */
  template<typename OtherDataT>
  MatrixXc operator*(const SpinOperator<OtherDataT>& other) const
  {
    const int aDim = X.rows();
    const int bDim = other.X.rows();
    const int dim = aDim * bDim;
    MatrixXc ret(dim, dim);
    const MatrixXc B_X = other.X;
    const MatrixXc B_Y = other.Y;
    const MatrixXc B_Z = other.Z;
    for(int i = 0; i < aDim; ++i) {
      for(int j = 0; j < aDim; ++j) {
//         cout << "block: (" << i << ',' << j << ") i.e. from " << i * bDim << " to " << j * bDim << " = " << endl;
//         cout << X(i, j) * B_X << endl << " + " << endl << Y(i, j) * B_Y << endl << " + " << endl << Z(i, j) * B_Z << endl;
//         cout << " = " << endl << X(i, j) * B_X + Y(i, j) * B_Y + Z(i, j) * B_Z << endl;
        ret.block(i * bDim, j * bDim, bDim, bDim) = X(i, j) * B_X + Y(i, j) * B_Y + Z(i, j) * B_Z;
//         cout << " = " << endl << ret.block(i * bDim, j * bDim, bDim, bDim) << endl;
      }
    }
    return ret;
  }

  DataT X;
  DataT Y;
  DataT Z;
};

MatrixXc directProductExpansion(const MatrixXc& A, const MatrixXc& B)
{
  const int aDim = A.rows();
  const int bDim = B.rows();
  const int dim = aDim * bDim;
  MatrixXc ret(dim, dim);
  for(int i = 0; i < aDim; ++i) {
    for(int j = 0; j < aDim; ++j) {
      ret.block(i * bDim, j * bDim, bDim, bDim) = A(i, j) * B;
    }
  }
  return ret;
}

/**
 * matrix product between 3x3 matrix and the Spin Operator (a tensor)
 */
template<typename DataT>
SpinOperator<DataT> operator*(const Matrix3c& A, const SpinOperator<DataT>& I)
{
    return SpinOperator<DataT>(
        A(0, 0) * I.X + A(0, 1) * I.Y + A(0, 2) * I.Z,
        A(1, 0) * I.X + A(1, 1) * I.Y + A(0, 2) * I.Z,
        A(2, 0) * I.X + A(2, 1) * I.Y + A(2, 2) * I.Z
    );
}

#undef ASSERT

#endif // MW_BACHELOR_OPERATORSUM_P_H
