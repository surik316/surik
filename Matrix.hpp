#ifndef MATRIX
#define MATRIX

#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>

template <typename T>
class Matrix
{
private:
    std::vector<std::vector<T>> matrix;
    size_t rowsCount = 0;
    size_t colsCount = 0;
public:
    Matrix<T>();
    Matrix<T>(size_t, size_t);
    T getEl(size_t, size_t);
    void changeEl(size_t, size_t, T);
    Matrix<T> transpose();
    Matrix<T> RowEchelonForm();
    Matrix<T> ReducedRowEchelonForm();
    Matrix<T> inverse();
    T det();
    size_t rank();
    size_t rows();
    size_t columns();
    const size_t columns() const;
    const size_t rows() const;
    std::vector<T>& operator[](size_t row) { return matrix[row]; }
    const std::vector<T>& operator[](size_t row) const { return matrix[row]; }
    void operator= (const Matrix<T>);

    Matrix& operator*=(const T& rhs) {
        for (auto& row : matrix) {
            for (auto& cell : row) {
                cell *= rhs;
            }
        }
        return *this;
    }

    Matrix& operator*=(const Matrix& rhs) {
	std::ofstream out("out.txt");
	if (colsCount != rhs.rowsCount) {
		if (out.is_open()) {
		  out << "-1" << std::endl;
		  out << "Cannot multiply matrixes with different sizes of rows and cols";
	          out.close();
		  exit(1);
		}
	}
        Matrix temp(rowsCount, rhs.colsCount);
        for (size_t i = 0; i < temp.rowsCount; i++) {
            for (size_t j = 0; j < temp.colsCount; j++) {
                for (size_t k = 0; k < colsCount; k++) {
                    temp[i][j] += matrix[i][k] * rhs[k][j];
                }
            }
        }
        std::swap(matrix, temp.matrix);
        std::swap(colsCount, temp.colsCount);
        return *this;
    }

    Matrix& operator+=(const Matrix& rhs) {
	std::ofstream out ("out.txt");
	if (rowsCount != rhs.rowsCount || colsCount != rhs.colsCount) {
		if (out.is_open()) {
		  out << "-1" << std::endl;
	          out << "Cannot add matrixes with different sizes of rols and cols";
		  out.close();
		  exit(1); 
		}
	}
        for (size_t i = 0; i < rowsCount; i++) {
            for (size_t j = 0; j < colsCount; j++) {
                matrix[i][j] += rhs[i][j];
            }
        }
    return *this;
    }
    
};

template <typename T>
Matrix<T>::Matrix() {}

template <typename T>
Matrix<T>::Matrix(size_t rows, size_t cols) {
    rowsCount = rows;
    colsCount = cols;
    for (unsigned int i = 0; i < rowsCount; i++) {
	    std::vector<T> tmp;
	    matrix.push_back(tmp);
	    for (unsigned int j = 0; j < colsCount; j++) {
		matrix[i].push_back(0);
	    }
    }
}
template <typename T>
void Matrix<T>::changeEl(size_t rCount, size_t cCount, T value) {
    if (rCount <= 0 || rCount > rowsCount) {
      exit(1);
    }
    if (cCount <= 0 || cCount > colsCount) {
      exit(1);
    }
    matrix[rCount-1][cCount-1] = value;
}

template <typename T>
T Matrix<T>::getEl(size_t rCount, size_t cCount) {
    if (rCount <= 0 || rCount > rowsCount) {
      exit(1);
    }
    if (cCount <= 0 || cCount > colsCount) {
      exit(1);
    }
    return matrix[rCount-1][cCount-1];
}


template <typename T>
T Matrix<T>::det() {
    unsigned int rowCount, colCount, step;
    T ret;
    Matrix A(rowsCount, colsCount);
    for (unsigned int i = 0; i < rowsCount; i++)
    for (unsigned int j = 0; j < colsCount; j++)
    A.matrix[i][j]=matrix[i][j];

    for(step = 0; step < rowsCount - 1; step++) {
        for(rowCount = step+1; rowCount < rowsCount; rowCount++) {
            ret = A.matrix[rowCount][step] / A.matrix[step][step];
            for (colCount = step; colCount < rowsCount; colCount++)
            A.matrix[rowCount][colCount] -= ret*A.matrix[step][colCount];
        }
    }
    ret = 1;
	for(rowCount = 0; rowCount < rowsCount; rowCount++) {
        ret = ret * A.matrix[rowCount][rowCount];
	}
    if (ret == -0) {ret = 0;}
    return ret;
}

template <typename T>
Matrix<T> Matrix<T>::transpose() {
    std::ofstream out ("out.txt");
    Matrix<T> ret(rowsCount, colsCount);
    if (rowsCount != colsCount) {
	if (out.is_open()) {
	   out << "-1" << std::endl;
	   out << "Cannot transpose non-square matrix";
	   out.close();
	   exit(1);
	}
    }
    for (unsigned int i = 0; i < rowsCount; i++) {
      	for (unsigned int j = 0; j < colsCount; j++) {
		ret.changeEl(i+1, j+1, matrix[j][i]);
	}
    }
    return ret;
}

template <typename T>
Matrix<T> Matrix<T>::RowEchelonForm() {
	Matrix<T> A = *this;
	unsigned int num = rowsCount > colsCount ? colsCount : rowsCount;
	T max, one;
	std::vector<T> temp;
	for (unsigned int i = 0; i < num; i++) {
		max = abs(A.matrix[i][i]);
		for (unsigned int j = i; j < rowsCount; j++) {
			if (abs(A.matrix[j][i]) > max) {
				max = abs(A.matrix[j][i]);
			}
		}
		if (max == 0) continue;

		for (unsigned int j = i + 1; j < rowsCount; j++) {
			one = A.matrix[j][i] / A.matrix[i][i];
			for (unsigned int k = 0; k < colsCount; k++) {
				A.matrix[j][k] -= A.matrix[i][k] * one;
			}
		}

	}
return A;
}

template <typename T>
Matrix<T> Matrix<T>::inverse() {
    int rows = rowsCount; //not unsigned
    int cols = colsCount; //not unsigned
    Matrix B(rows,cols),  C(cols,cols*2); // with E matrix
    int row, column, step;
    T ret;
    std::ofstream out ("out.txt");

    if(rows != cols) {
	if (out.is_open()) {
		out << "-1" << std::endl;
		out << "Cannot inverse non-square matrix"; 
		out.close();
		exit(1);
	}
    }

    for(row = 0; row < rows; row++)
    for(column = 0; column < cols; column++)
    B.matrix[row][column] = matrix[row][column]; 

    if (B.det() == 0) {
	    if (out.is_open()) {
		    out << "-1" << std::endl;
		    out << "Cannot inverse matrix with det = 0";
		    out.close();
		    exit(1);
	    }
    }

    for(row = 0; row < rows; row++)
    for(column = 0; column < rows; column++)
    C.matrix[row][column] = matrix[row][column];

    for(int i = 0; i < rows; i++)
    C.matrix[i][rows+i] = 1;

    for(step = 0; step < rows - 1; step++) {
        for(row = step + 1; row < rows; row++) {
            ret = C.matrix[row][step] / C.matrix[step][step];
            for (column = step; column < 2 * rows; column++)
            C.matrix[row][column] -= ret*C.matrix[step][column];
      }
    }
    for(step = 1; step <= rows - 1; step++) {
        for(row = rows - step -1; row >= 0; row--) {
            ret = C.matrix[row][rows-step] / C.matrix[rows-step][rows-step];
            for (column = rows; column < 2 * rows; column++)
            C.matrix[row][column] -= ret*C.matrix[rows-step][column];
        }
    }

    for(row = 0; row < rows; row++)
    for(column = 0; column < rows; column++)
    B.matrix[row][column] = C.matrix[row][rows+column] / C.matrix[row][row];
return B;
}

template <typename T>
Matrix<T> Matrix<T>::ReducedRowEchelonForm() {
    Matrix<T> A = *this;
    int iom;
    std::vector<T> t;
    A.RowEchelonForm();
    T one, second;
    for (unsigned int i = 0; i < rowsCount; i++) {
        iom = -1;
        for (unsigned int j = 0; j < colsCount; j++) {
            if (A.matrix[i][j] != 0) {
                iom = j;
                break;
            }
        }
        if (iom == -1) continue;

        one = 1/A.matrix[i][iom];

        for (unsigned int j = 0; j < colsCount; j++) {
            A.matrix[i][j] *= one;
            if (A.matrix[i][j] == -0) A.matrix[i][j] = 0;
        }

        for (unsigned int j = 0; j < rowsCount; j++) {
            if (i == j) continue;
            second = A.matrix[j][iom];
            for (unsigned int k = 0; k < colsCount; k++) {
                A.matrix[j][k] -= A.matrix[i][k]*second;
            }
        }
    }
    for (unsigned int j = 0; j < rowsCount - 1; j++) {
            for (unsigned int k = 0; k < rowsCount - j - 1; k++) {
                int it1 = colsCount, it2 = colsCount;
                for (unsigned int z = 0; z < colsCount; z++) {
                    if (A.matrix[k][z] != 0) {
                        it1 = z;
                        break;
                    }
                }
                for (unsigned int z = 0; z < colsCount; z++) {
                    if (A.matrix[k + 1][z] != 0) {
                        it2 = z;
                        break;
                    }
                }
                if (it2 < it1) {
                   t = A.matrix[k + 1];
                   A.matrix[k + 1] = A.matrix[k];
                   A.matrix[k] = t;
                }
            }
        }

    return A;
}

template <typename T>
size_t Matrix<T>::rank() {
    size_t ret = 0;
    Matrix<T> A = *this;
    A = A.ReducedRowEchelonForm();
    for (unsigned int i = 0; i < A.rows(); i++) {
       for (unsigned int j = 0; j < A.columns(); j++) {
          if (A.getEl(i+1, j+1) != 0) {
		  ret++;
		  break;
	  }
       }
    }
return ret; 
}

template <typename T>
Matrix<T> operator+ (Matrix<T> lhs, const Matrix<T>& rhs) { return lhs += rhs; }

template <typename T>
Matrix<T> operator* (Matrix<T> lhs, const Matrix<T>& rhs) { return lhs *= rhs; }

template <typename T>
Matrix<T> operator* (Matrix<T> lhs, const T& rhs) { return lhs *= rhs; }

template <typename T>
Matrix<T> operator* (const T& lhs, Matrix<T> rhs) { return rhs *= lhs; }

template <typename T>
size_t Matrix<T>::rows() { return rowsCount; }

template <typename T>
size_t Matrix<T>::columns() {return colsCount; }

template <typename T>
const size_t Matrix<T>::columns() const { return colsCount; }

template <typename T>
const size_t Matrix<T>::rows() const { return rowsCount; }

template <typename T>
void Matrix<T>::operator= (const Matrix<T> rhs) {
    rowsCount = rhs.rowsCount;
    colsCount = rhs.colsCount;
    matrix = rhs.matrix;
}

template <typename T>
std::istream& operator >> (std::istream& is, Matrix<T>& mat) {
    for (size_t i = 0; i < mat.rows(); i++) {
        for (size_t j = 0; j < mat.columns(); j++) {
            is >> mat[i][j];
        }
    }
return is;
}

template <typename T>
std::ostream& operator<< (std::ostream& os, const Matrix<T>& matrix) {
    for (size_t i = 0; i < matrix.rows(); i++) {
        for (size_t j = 0; j < matrix.columns(); j++) {
            os << matrix[i][j] << ' ';
        }
        os << "\n";
    }
return os;
}

#endif
