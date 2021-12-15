/*
 *------------------------------------------------------
 * (c) 2021, mohab metwally All Rights Reserved
 *
 * the content of this file may not be disclosed, copied, or duplicated in
 * any form, in whole, or part, without the prior wirtten permission of
 *  mohab metwally
 *------------------------------------------------------
 *
 * @author: mohab metwally
 * @date: 2021/12/1
*/

#pragma once

#include <iostream>
#include <memory>
#include <random>

template<class T, unsigned int N>
class Point;

enum MatrixInit 
{
  kIdentity = 0,
  kZeros,
  kOnes,
  kRandom
};

/**
  * class Matrix
  * This class is a general Matrix type that is also a base for special matrices like Transformation
  */
template<class T>
class Matrix
{
public:

  // Constructors/Destructors
  //  

  /**
  * This constructor initializes the Matrix with one of the types Zeros, Ones, Identity, Random
  * @param  rows: number of rows
  * @param  cols: number of columns
  * @param  init_type: initial value (kZeros, kOnes, kIdentity, kRandom
  */
  Matrix(const int rows = 1, const int cols = 1, const MatrixInit init_type = kZeros) : rows_(rows), cols_(cols)
  {
    if (init_type == kRandom)
    {
      matrix_.resize(rows_, std::vector<T>(cols, 0));

      std::random_device rd;
      std::mt19937 mt(rd());
      std::uniform_real_distribution<T> dist(0, 1);
      for (int r = 0; r < rows_; ++r) {
        for (int c = 0; c < cols_; ++c) {
          matrix_[r][c] = dist(mt);
        }
      }
    }
    else if (init_type == kZeros) {
      matrix_.resize(rows_, std::vector<T>(cols, 0));
    }
    else if (init_type == kOnes) {
      matrix_.resize(rows_, std::vector<T>(cols, 1));
    }
    else if (init_type == kIdentity) {
      matrix_.resize(rows_, std::vector<T>(cols, 0));
      for (int i = 0; i < std::min(rows_, cols_); ++i) {
        matrix_[i][i] = 1;
      }
    }
  }

  /**
  * This is the Copy Constructor
  */
  Matrix(const Matrix& matrix_b) : rows_(matrix_b.rows_), cols_(matrix_b.cols_)
  {
    matrix_ = matrix_b.matrix_;
  }

  /**
  * This is the assignment operator
  * Copy the Matrix at the rhs of the oeprator = to the current matrix
  * @param matrix_b: the matrix at the right hand side of the = operator
  * returns reference to the current object (used for nested = operations)
  */
  Matrix<T>& operator=(const Matrix<T>& matrix_b) {
    rows_ = matrix_b.rows_;
    cols_ = matrix_b.cols_;
    matrix_ = matrix_b.matrix_;
    return *this;
  }

  /**
  * This is a constructor with 1-D vector in row-major order
  */
  Matrix(const int rows, const int cols, const std::vector<T>& data) : rows_(rows), cols_(cols)
  {
    if (data.size() != rows * cols) {
      throw "data size must equal rows x cols";
    }
    for (int r = 0; r < rows_; ++r) {
      matrix_.push_back(std::vector<T>(data.begin() + r * cols_, data.begin() + (r + 1) * cols_));
    }

  }

  /**
* This is a constructor with 1-D array in row-major order
*/
  Matrix(const int rows, const int cols, T* data) : rows_(rows), cols_(cols)
  {
    
    for (int r = 0; r < rows_; ++r) {
      matrix_.push_back(std::vector<T>(data.begin() + r * cols_, data.begin() + (r + 1) * cols_));
    }
  }


  /**
  * Empty Destructor
  */
  virtual ~Matrix()
  {
  }


  /**
  * This fucntion gets the matrix data in 1-D vector in row-major order
  */
  std::vector<T> GetData()
  {
    std::vector<T> res;

    for (int r = 0; r < rows_; ++r) {
      res.insert(res.end(), matrix_[r].begin(), matrix_[r].end());
    }
    return res;
  }

  bool operator==(const Matrix& matrix_b) const
  {
    for (unsigned int r = 0; r < rows_; ++r) {
      for (unsigned int c = 0; c < cols_; ++c) {
        if (matrix_[r][c] != matrix_b.matrix_[r][c]) {
          return false;
        }
      }
    }
    return true;
  }


  /**
   * This function multiplies matrix with another matrix
   * @return Matrix
   * @param  matrix_b
   */
  Matrix operator*(const Matrix& matrix_b) const
  {
    if (cols_ != matrix_b.rows_) {
      throw "cols(matrix_a) must equal rows(matrix_b) to multiply!";
    }
    Matrix<T> res(rows_, matrix_b.cols_);

    for (unsigned int k = 0; k < matrix_b.cols_; ++k) {
      for (unsigned int i = 0; i < rows_; ++i) {
        T tmp = 0;
        for (unsigned int j = 0; j < cols_; ++j) {
          tmp += matrix_[i][j] * matrix_b.matrix_[j][k];
        }
        res[i][k] = tmp;
      }
    }
    return res;
  }

  /**
   * This function multiplies the current matrix with vector/point
   * @return Matrix
   * @param  point
   */
  template<unsigned int N>
  Point<T, N> operator*(Point<T, N> point) const
  {
    if (N != GetColumns()) {
      throw "cols(matrix_a) must equal point size to multiply!";
    }
    Point<T, N> res;

    for (unsigned int i = 0; i < rows_; ++i) {
      T tmp = 0;
      for (unsigned int j = 0; j < cols_; ++j) {
        tmp += matrix_[i][j] * point[j];
      }
      res[i] = tmp;
    }
    return res;
  }


  /**
   * This function multiplies the current matrix with a scalar value
   * @return Matrix
   * @param  number
   */
  //template<class type>
  Matrix operator*(const T number)
  {
    Matrix<T> res(rows_, cols_);
    for (unsigned int r = 0; r < rows_; ++r) {
      for (unsigned int c = 0; c < cols_; ++c) {
        res.matrix_[r][c] = matrix_[r][c] * number;
      }
    }
    return res;
  }


  /**
   * This function divides the current matrix with a scalar
   * @return Matrix
   * @param  denominator
   */
  template<class type>
  Matrix operator/(const type denominator) const
  {
    if (denominator == 0) {
      throw "Divison by Zero!";
    }
    Matrix<T> res(rows_, cols_);
    for (unsigned int r = 0; r < rows_; ++r) {
      for (unsigned int c = 0; c < cols_; ++c) {
        res[r][c] = matrix_[r][c] / denominator;
      }
    }
    return res;
  }

  /**
 * This function divides the current matrix with a scalar
 * @param  denominator
 */
  template<class type>
  void operator/=(const type denominator)
  {
    if (denominator == 0) {
      throw "Divison by Zero!";
    }
    for (unsigned int r = 0; r < rows_; ++r) {
      for (unsigned int c = 0; c < cols_; ++c) {
        matrix_[r][c] /= denominator;
      }
    }
  }

  /**
   * This function finds the the power of a matrix
   * @return Matrix
   * @param  exponent
   */
  Matrix operator^(int exponent) const
  {
    if (rows_ != cols_) {
      throw "power is defined only for square Matrices";
    }
    Matrix res(rows_, cols_, kIdentity);
    Matrix tmp(*this);
    if (exponent < 0) {
      exponent = -exponent;
      tmp = tmp.GetInverse();
    }

    while (exponent) {
      if (exponent & 1) {
        res = res * tmp;
      }
      exponent >>= 1;
      tmp = tmp * tmp;
    }
    return res;
  }


  /**
   * This function adds two matrices
   * @return Matrix
   * @param  matrix_b
   */
  Matrix operator+(const Matrix& matrix_b) const
  {
    if (rows_ != matrix_b.rows_ || cols_ != matrix_b.cols_) {
      throw "Matrices must be of same size to Add";
    }
    Matrix<T> res(rows_, cols_);
    for (unsigned int r = 0; r < rows_; ++r) {
      for (unsigned int c = 0; c < cols_; ++c) {
        res[r][c] = matrix_[r][c] + matrix_b.matrix_[r][c];
      }
    }
    return res;
  }

  /**
 * This function adds a matrix to the current matrix
 * @return Matrix
 */
  void operator+=(const Matrix & matrix_b)
  {
    if (rows_ != matrix_b.rows_ || cols_ != matrix_b.cols_) {
      throw "Matrices must be of same size to Add";
    }
    for (unsigned int r = 0; r < rows_; ++r) {
      for (unsigned int c = 0; c < cols_; ++c) {
        matrix_[r][c] += matrix_b.matrix_[r][c];
      }
    }
  }


  /**
   * This function Subtracts two matrices
   * @return Matrix
   * @param  matrix_b
   */
  Matrix operator-(const Matrix& matrix_b) const
  {
    if (rows_ != matrix_b.rows_ || cols_ != matrix_b.cols_) {
      throw "Matrices must be of same size to Subtract";
    }
    Matrix<T> res(rows_, cols_);
    for (unsigned int r = 0; r < rows_; ++r) {
      for (unsigned int c = 0; c < cols_; ++c) {
        res[r][c] = matrix_[r][c] - matrix_b.matrix_[r][c];
      }
    }
    return res;
  }

  /**
   * This function Negates the current matrix an return a copy of the negated matrix
   * @return Matrix
   */
  Matrix operator-()
  {
    Matrix res = *this * -1;

    return res;
  }


  /**
  * This function sets one value inside the matrix by row and column
  * @return true if the set operation was successful
  * @param  row
  * @param  column
  * @param  value
  */
  inline bool SetElement(const int row, const int column, const T value)
  {
    if (row >= 0 && row < rows_ && column >= 0 && column < cols_) {
      matrix_[row][column] = value;
      return true;
    }
    return false;
  }

  /*
  * This function finds and returns the minimum coefficient in the Matrix
  */
  T GetMin() const
  {
    T res = matrix_[0][0];
    for (unsigned int r = 0; r < rows_; ++r) {
      for (unsigned int c = 0; c < cols_; ++c) {
        if (res > matrix_[r][c]) {
          res = matrix_[r][c];
        }
      }
    }

    return res;
  }

  /*
  * This function finds and returns the maximum coefficient in the Matrix
  */
  T GetMax() const
  {
    T res = matrix_[0][0];
    for (unsigned int r = 0; r < rows_; ++r) {
      for (unsigned int c = 0; c < cols_; ++c) {
        if (res < matrix_[r][c]) {
          res = matrix_[r][c];
        }
      }
    }

    return res;
  }

  /*
  * This function estiamtes and returns the mean of the coefficient in the Matrix
  */
  T GetMean() const
  {
    T res = 0;
    for (unsigned int r = 0; r < rows_; ++r) {
      for (unsigned int c = 0; c < cols_; ++c) {
        res += matrix_[r][c];
      }
    }

    return res / (rows_ * cols_);
  }

  /*
  * This function estimates and returns the trace of the Matrix (in case it is square)
  */
  T GetTrace() const
  {
    if (rows_ != cols_) {
      throw "Trace is defined only for square Matrix";
    }
    T result = 0;
    for (unsigned int i = 0; i < rows_; ++i) {
      result += matrix_[i][i];
    }

    return result;
  }

  /**
  * This function extracts and returns a sub matrix identified by the given parameters
  * @return Matrix
  * @param  start_row
  * @param  start_column
  * @param  rows
  * @param  columns
  */
  Matrix GetSubMatrix(const int start_row, const int start_column, const int num_rows, const int num_cols) const
  {
    Matrix<T> res(num_rows, num_cols);
    for (unsigned int r = 0; r < num_rows; r++) {
      for (unsigned int c = 0; c < num_cols; c++) {
        res.matrix_[r][c] = matrix_[(start_row + r)][(start_column + c)];
      }
    }
    return res;
  }

  /**
  * This functionestimates the inverse of the current matrix and returns this inverse
  * @return Matrix inverse of the current matrix if exists
  */
  Matrix GetInverse() const
  {
    if (rows_ != cols_) {
      throw "Inverse needs square matrix";
    }
    if (Determinant() == 0) {
      throw "Can not find Inverse for this Matrix, Determinant = 0!";
    }
    Matrix id(rows_, cols_, kIdentity);

    auto matrix(*this);

    for (int i = 0; i < rows_; ++i) {
      // try to eliminate all column nunmber i except matrix[i][i]
      T pivot = matrix[i][i];
      if (pivot == 0) {
        // find another pivot
        for (unsigned int j = i + 1; j < rows_; ++j) {
          if (matrix[j][i] != 0) {
            pivot = matrix[j][i];
            std::swap(matrix[i], matrix[j]);
            std::swap(id[i], id[j]);
            break;
          }
        }
      }
      for (unsigned int r = 0; r < rows_; ++r) {
        if (r == i) continue;
        T tmp = matrix[r][i];

        if (tmp != 0) {
          double ratio = tmp / pivot;
          // now eliminate
          for (int c = 0; c < cols_; ++c) {
            matrix[r][c] -= ratio * matrix[i][c];
            // change the identity as well
            id[r][c] -= ratio * id[i][c];
          }
        }
      }
    }
    // divide by the values in the diagonal
    for (unsigned int r = 0; r < rows_; ++r) {
      for (unsigned int c = 0; c < cols_; ++c) {
        id[r][c] /= matrix[r][r];
      }
    }
    return id;

  }


  /**
  * This function estimates the transpose of the current matrix in a new matrix
  * @return Matrix transpose of the current matrix
  */
  Matrix Transpose() const
  {
    Matrix res(cols_, rows_);
    for (unsigned int r = 0; r < rows_; ++r) {
      for (unsigned int c = 0; c < cols_; ++c) {
        res[c][r] = matrix_[r][c];
      }
    }
    return res;
  }

  /**
  * @return int number of rows
  */
  int GetRows() const
  {
    return rows_;
  }

  /**
  * @return int number of columns
  */
  int GetColumns() const
  {
    return cols_;
  }

  /**
  * This function gets the element specified by the row and column
  * @return T
  * @param  row
  * @param  column
  */
  inline T GetElement(int row, int column) const
  {
    if (row < 0 || row >= rows_ || column < 0 || column >= cols_) {
      throw "Matrix index out of range!";
    }
    return matrix_[row][column];
  }

  /**
  * This function finds and returns the determinant of the current matrix if it was a square
  * @return T the determinant of the current matrix
  */
  T Determinant() const
  {
    if (rows_ != cols_) {
      throw "Matrix must be square to find Determinant!";
    }
    if (rows_ == 1) {
      return matrix_[0][0];
    }
    auto matrix(*this);
    int swaps = 0;
    for (unsigned int i = 0; i < rows_; ++i) {
      // try to eliminate all column nunmber i except matrix[i][i]
      T pivot = matrix[i][i];
      if (pivot == 0) {
        // find another pivot
        for (unsigned int j = i + 1; j < rows_; ++j) {
          if (matrix[j][i] != 0) {
            pivot = matrix[j][i];
            std::swap(matrix[i], matrix[j]);
            ++swaps;
            break;
          }
        }
        if (pivot == 0) {
          return 0;
        }
      }
      for (unsigned int r = i + 1; r < rows_; ++r) {
        T tmp = matrix[r][i];

        if (tmp != 0) {
          double ratio = tmp / pivot;
          // now eliminate
          for (unsigned int c = 0; c < cols_; ++c) {
            matrix[r][c] -= ratio * matrix[i][c];
          }
        }
      }
    }
    // determinant now is the multiplication of the diagonal
    T res = matrix[0][0];
    for (unsigned int i = 1; i < rows_; ++i) {
      res *= matrix[i][i];
    }
    // check the sign based on how many swaps
    if (swaps & 1) {
      return -res;
    }
    return res;
  }

  /**
  * This function returns the row spefified by row_index by reference
   * @return std::vector<T>& the row specified by the given row_index
   * @param  row_inedx the row index to be returned
   */
  std::vector<T>& operator[](int row_index) {
    if (row_index >= 0 && row_index < rows_) {
      return matrix_[row_index];
    }
    throw "row_index out of matrix indeces";
  }


  /* 
  * Sets matrix column values
  * @param colValues: column values
  * @param num: column number
  */ 
  void SetCol(const std::vector<T>& colValues, unsigned int num) 
  {
    if(colValues.size() == rows_ || num < cols_)
    {
      for (int r = 0; r < rows_; ++r) {
        SetElement(r, num, colValues[r]);
      }
    }
  }

protected:
  std::vector<std::vector<T> > matrix_; // the matrix data represented by vector of vectors in row major order
  unsigned int rows_; // number of rows in the current matrix
  unsigned int cols_; // number of columns in the current matrix

};


/**
  * class Transformation
  * This class is a special class inherits the matrix class with matrix size 4x4
  */

template<class T>
class Transformation : public Matrix<T>
{
public:

  /**
   * Empty Constructor
   */ 
  Transformation() : Matrix<T>(4, 4, kIdentity)
  {
  }

  /*
  * This is the copy constructor
  */
  Transformation(Transformation& matrix) : Matrix<T>(4, 4) 
  {
    this->matrix_ = matrix.matrix_;
  }

  /**
 * This is a constructor with 1-D vector in row-major order
 */
  Transformation(std::vector<T> &data) : Matrix<T>(4, 4, data)
  {

  }

  /*
  * Transformation constructor from the base class Matrix
  */
  Transformation(Matrix<T> &matrix) : Matrix<T>(4, 4, kIdentity)
  {
	  if (matrix.GetColumns() != 4 || matrix.GetRows() != 4) {
		  throw "matrix expected to be 4x4 transformation matrix";
	  }
	  Transformation(matrix.GetData());
  }

  /**
   * Empty Destructor
   */
  virtual ~Transformation()
  {
  }
private:
};
