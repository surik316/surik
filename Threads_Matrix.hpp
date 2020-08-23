//
//  Matrix.hpp
//  
//
//  Created by Максим Сурков on 31.05.2020.
//
#ifndef MATRIX
#define MATRIX
#include <iostream>
#include <fstream>
#include <thread>
#include <mutex>
#include <vector>
#include <future>
#include <exception>
using namespace std;
template<typename T>
class Matrix{
private:
    static size_t parallel;
    size_t rows;
    size_t columns;
    vector<vector<T> > matrix;
    
public:
    
    T Multiply(const vector <T>&, const vector <T>& );
    Matrix<T> Transpose() const;
    Matrix() : rows(0), columns(0){}
    void LoadMatrix(const char* );
    void LoadMatrix(ifstream& );
    void SaveMatrix(const char *) const;
    
    Matrix(size_t r, size_t c) : rows(r), columns(c){
        matrix.resize(rows);
        for (size_t i = 0;  i < rows; ++i){
            matrix[i].resize(columns);
        }
    }
    
    static void SetParallel(size_t &num){
        parallel = num;
    }
    
    static size_t GetParallel(){
        return parallel;
    }
    
    Matrix operator * ( const Matrix& rhs){
        
        if (columns != rhs.rows){
            throw std::invalid_argument("Can't multiply");
        }
        
        Matrix<T> temp(rows, rhs.columns);
        size_t counter = 0;
        vector<future<int>> futures(parallel);
        vector<promise<int>> promises(parallel);
        Matrix F = rhs.Transpose();
        mutex door;
        
        auto lambda  = [&](size_t x){
            
            while (true){
                door.lock();
                size_t cnt = counter;
                counter++;
                door.unlock();
                
                if (cnt >= rows * rhs.columns){
                    break;
                }
                
                temp.matrix[cnt / rhs.columns][cnt % rhs.columns ] = Multiply(matrix[cnt / rhs.columns], F.matrix[cnt % rhs.columns]);
            }
            promises[x].set_value_at_thread_exit(1);
        };
        
        for (size_t i = 0; i < futures.size() ; ++i){
            futures[i] = promises[i].get_future();
        }
        
        for (size_t i =0 ; i < futures.size(); ++i){
            thread(lambda, i).detach();
        }
        
        for (size_t i =0 ; i < futures.size(); ++i){
            futures[i].wait();
        }
        
        for (size_t i =0 ; i < futures.size(); ++i){
             cout << "thread "<< i+1 <<" completed his work" << endl;
        }
        
        return temp;
    }
};

template<typename T>
size_t Matrix<T>::parallel = 0;

template<typename T>
Matrix<T> Matrix<T>::Transpose() const{
    Matrix temp(columns, rows);
    
    for (size_t i =0; i < rows; ++i){
        
        for (size_t j =0; j < columns; ++j){
            temp.matrix[j][i] = matrix[i][j];
        }
        
    }
    return temp;
}

template<typename T>
void Matrix<T>::LoadMatrix(const char* filename){
    ifstream fin(filename);
    LoadMatrix(fin);
}

template<typename T>
void Matrix<T>::LoadMatrix(ifstream& fin){
    size_t r, c;
    fin >> r >> c;
    rows = r;
    columns = c;
    
    matrix.resize(rows);
    for (size_t i = 0 ; i < rows; ++i){
        matrix[i].resize(columns);
        for (size_t j = 0; j < columns; ++j){
            fin >> matrix[i][j];
        }
    }
}

template<typename T>
void Matrix<T>::SaveMatrix(const char * filename) const{
    ofstream fout(filename);
    fout << rows << " " << columns << endl;
    
    for (size_t i = 0; i < rows; ++i) {
        
        for (size_t j = 0; j < columns; ++j) {
            fout << matrix[i][j] << " ";
        }
        
        fout << endl;
    }
    fout.close();
}

template<typename T>
T Matrix<T>::Multiply(const vector <T>& v1,const  vector <T>& v2){
    T ret = 0;
    
    for (size_t i = 0; i < v1.size(); ++i){
        ret += v1[i] * v2[i];
    }
    
    return ret;
}
#endif
