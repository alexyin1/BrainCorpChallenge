#pragma once
#include <thread>
#include <vector> 
#include <string>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>       /* fabs */

namespace MatOps  // namespace declaration
{
    // Global vars
    double EPSILON = 1e-5;
    std::string TESTS_DIR = "tests/";
    int N_MM_TESTS = 15;
    int N_TP_TESTS = 11;

    template <class Dtype>
    class Matrix
    {
        int nRows;
        int nCols;
        std::vector<Dtype> data;

    public:
        //empty data constructor
        Matrix(int nRows, int nCols){
            this->nRows = nRows;
            this->nCols = nCols;
        }
        //copy constructor
        Matrix(const Matrix<Dtype> &m){
            this->nRows = m.nRows;
            this->nCols = m.nCols;
            this->data = m.data;
        }
        // vec constructor
        Matrix(std::vector<Dtype> arr, int nRows, int nCols){
            this->nRows = nRows;
            this->nCols = nCols;
            this->data = arr;
        }
        // 1D array constructor
        Matrix(Dtype* arr, int nRows, int nCols){
            this->nRows = nRows;
            this->nCols = nCols;
            for (int i=0; i<nRows*nCols; i++){
                this->data.push_back(arr[i]);
            }
        }
        // 2D array constructor
        Matrix(Dtype** arr, int nRows, int nCols){
            this->nRows = nRows;
            this->nCols = nCols;
            for (int i=0; i<nRows; i++){
                for (int j=0; j<nCols; j++){
                    this->data.push_back(arr[i][j]);
                }
            }
        }

        int getRows() const{
            return nRows;
        }

        int getCols() const{
            return nCols;
        }
        std::vector<Dtype> getData() const{
            return this->data;
        }

        Dtype at(int row, int col) const{
            return this->data[row*this->nCols+col];
        }

        // TODO: add try catch to handle type(A) == type(B), ensure
        // TODO: multiprocess matmul via blocks
        //returns C, C=A*B;
        Matrix operator*(const Matrix& B){
            assert (this->nCols==B.nRows);
            int nRows = this->nRows;
            int nCols = B.nCols;
            int aCols = this->nCols;
            std::vector<Dtype> C;
            C.reserve(nRows*nCols);
            int n_dot = this->nCols;
            for (int i=0; i<nRows; i++){
                for (int j=0; j<nCols; j++){
                    for (int k=0; k<n_dot; k++){
                        if (i*nCols+j < C.size()){
                            C[i*nCols+j] += this->data[i*aCols+k] * B.data[k*nCols+j];
                        }
                        else{
                            C.push_back(this->data[i*aCols+k] * B.data[k*B.nCols+j]);
                        }
                    }
                }
            }
            return Matrix(C, nRows, nCols);
        }

        Matrix matmul(const Matrix& b){
            return (*this)*b;
        }
        //TODO: in-place transpose of square matrices
        //TODO : in-place transpose of non-square matrices
        //void transposeInPlace(){

        //out of place transpose
        Matrix transpose(){
            std::vector<Dtype> temp;
            temp.reserve(this->nRows*this->nCols);
            for (int i=0; i<this->nCols; i++){
                for (int j=0; j<this->nRows; j++){
                    temp.push_back(this->at(j,i));
                }
            }
            return Matrix<Dtype>(temp, this->nCols, this->nRows);
        };
    };

    template <class Dtype>
    void print(Matrix<Dtype> a){
        std::cout << a;        
    };

    template <class Dtype>
    std::ostream& operator<<(std::ostream &out, const Matrix<Dtype> &m)
    {
        int nRows = m.getRows();
        int nCols = m.getCols();
        std::vector<Dtype> data = m.getData();
        out << typeid(Dtype).name() << "," << nRows << ",";
        out << nCols << "\n";
        for (int i =0; i <nRows; i++){
            for (int j=0; j<nCols; j++){
                out << data[nCols*i+j] << ",";
                if (j == nCols-1)
                    out << "\n";
            }
        }
        return out;
    }

    void readIntMat(std::fstream& f, std::vector<int>& v){
        std::string line;
        if (f.is_open()){
            while (std::getline(f,line)){
                if (line.empty()){
                    break;
                }
                else{
                    std::stringstream ss(line);
                    while (ss.good() && ss.peek()!='\n'){
                        std::string temp;
                        std::getline(ss, temp, ',');
                        if (!temp.empty()){
                            int val = std::stoi(temp);
                            v.push_back(val);
                        }
                    }
                }
            }
        }
        return;
    }

    void readFloatMat(std::fstream& f, std::vector<float>& v){
        std::string line;
        if (f.is_open()){
            while (std::getline(f,line)){
                if (line.empty()){
                    break;
                }
                else{
                    std::stringstream ss(line);
                    while (ss.good() && ss.peek()!='\n'){
                        std::string temp;
                        std::getline(ss, temp, ',');
                        if (!temp.empty()){
                            float val = std::stof(temp);
                            v.push_back(val);
                        }
                    }
                }
            }
        }
        return;
    }

    void readDoubleMat(std::fstream& f, std::vector<double>& v){
        std::string line;
        if (f.is_open()){
            while (std::getline(f,line)){
                if (line.empty()){
                    break;
                }
                else{
                    std::stringstream ss(line);
                    while (ss.good() && ss.peek()!='\n'){
                        std::string temp;
                        std::getline(ss, temp, ',');
                        if (!temp.empty()){
                            double val = std::stod(temp);
                            v.push_back(val);
                        }
                    }
                }
            }
        }
        return;
    }

    // account for floating point precision loss
    template <class Dtype>
    bool approx(Dtype a, Dtype b){
        return fabs(a - b) < EPSILON;
    }

    template <class Dtype>
    bool equals(Matrix<Dtype>& calc, Matrix<Dtype>& expected){
        int nRows = expected.getRows();
        int nCols = expected.getCols();
        if (calc.getRows() != nRows || calc.getCols() != nCols){
            return false;
        }
        for (int i =0; i <nRows; i++){
            for (int j=0; j<nCols; j++){
                if (!approx(calc.at(i,j), expected.at(i,j))) {
                    return false;
                }
            }
        }
        return true;
    }

    void readParams(std::fstream& f, int& rows, int& cols, std::string& dtype){
        std::string line;
        std::getline(f,line);
        std::string temp;
        std::stringstream ss(line);
        std::getline(ss, temp, ',');
        rows = std::stoi(temp);
        std::getline(ss, temp, ',');
        cols = std::stoi(temp);
        std::getline(ss, temp, ',');
        dtype = temp;
    }

    template<class Dtype>
    void printResult(Matrix<Dtype>& calc, Matrix<Dtype>& expected, std::string& name){
        if (equals(calc, expected)){
            printf("Test %s Passed!\n", name.c_str());
        }
        else{
            printf("Test %s Failed\n", name.c_str());
        }
    }

    void run_tpcase(std::string name){
        int rows;
        int cols;
        std::string dtype;
        std::fstream f;
        std::string fname = TESTS_DIR + name;
        f.open(fname);
        readParams(f, rows, cols, dtype);
        if (dtype == "int"){
            std::vector<int> a;
            std::vector<int> At;
            //read A;
            readIntMat(f, a);
            Matrix<int> A(a, rows, cols);
            //read At
            readParams(f, rows, cols, dtype);
            readIntMat(f, At);
            Matrix<int> expected(At, rows, cols);
            Matrix<int> calc(rows, cols);
            calc = A.transpose();
            printResult(calc, expected, name);
        }
        else if (dtype == "float"){
            std::vector<float> a;
            std::vector<float> At;
            //read A;
            readFloatMat(f, a);
            Matrix<float> A(a, rows, cols);
            //read At
            readParams(f, rows, cols, dtype);
            readFloatMat(f, At);
            Matrix<float> expected(At, rows, cols);
            Matrix<float> calc(rows, cols);
            calc = A.transpose();
            printResult(calc, expected, name);
        }
        else if (dtype == "double"){
            std::vector<double> a;
            std::vector<double> At;
            //read A;
            readDoubleMat(f, a);
            Matrix<double> A(a, rows, cols);
            //read At
            readParams(f, rows, cols, dtype);
            readDoubleMat(f, At);
            Matrix<double> expected(At, rows, cols);
            Matrix<double> calc(rows, cols);
            calc = A.transpose();
            printResult(calc, expected, name);
        }
        f.close();

    };


    void run_mmcase(std::string name){
        int rows;
        int cols;
        std::string dtype;
        std::fstream f;
        std::string fname = TESTS_DIR + name;
        f.open(fname);
        readParams(f, rows, cols, dtype);
        if (dtype == "int"){
            std::vector<int> a;
            std::vector<int> b;
            std::vector<int> c;
            //read A and B;
            readIntMat(f, a);
            Matrix<int> A(a, rows, cols);
            readParams(f, rows, cols, dtype);
            readIntMat(f, b);
            Matrix<int> B(b, rows, cols);
            //read C
            readParams(f, rows, cols, dtype);
            readIntMat(f, c);
            Matrix<int> expected(c, rows, cols);
            Matrix<int> calc(rows, cols);
            calc = A*B;
            printResult(calc, expected, name);
        }
        else if (dtype == "float"){
            std::vector<float> a;
            std::vector<float> b;
            std::vector<float> c;
            //read A and B;
            readFloatMat(f, a);
            Matrix<float> A(a, rows, cols);
            readParams(f, rows, cols, dtype);
            readFloatMat(f, b);
            Matrix<float> B(b, rows, cols);
            //read C
            readParams(f, rows, cols, dtype);
            readFloatMat(f, c);
            Matrix<float> expected(c, rows, cols);
            Matrix<float> calc(rows, cols);
            calc = A*B;
            printResult(calc, expected, name);
        }
        else if (dtype == "double"){
            std::vector<double> a;
            std::vector<double> b;
            std::vector<double> c;
            //read A and B;
            readDoubleMat(f, a);
            Matrix<double> A(a, rows, cols);
            readParams(f, rows, cols, dtype);
            readDoubleMat(f, b);
            Matrix<double> B(b, rows, cols);
            readParams(f, rows, cols, dtype);
            //read C
            readDoubleMat(f, c);
            Matrix<double> expected(c, rows, cols);
            Matrix<double> calc(rows, cols);
            calc = A*B;
            printResult(calc, expected, name);
        }
        f.close();
    };

    void runTests(){
        // matmul tests
        for (int i=0; i<N_MM_TESTS; i++){
            std::string prefix = "mm_";
            if (i < 10){
                prefix += "00" + std::to_string(i) + ".txt"; 
            }
            else if (i >= 10 && i < 100){
                prefix += "0" + std::to_string(i) + ".txt";
            }
            run_mmcase(prefix);
        }
        // transpose tests
        for (int i=0; i<N_TP_TESTS; i++){
            std::string prefix = "tp_";
            if (i < 10){
                prefix += "00" + std::to_string(i) + ".txt"; 
            }
            else if (i >= 10 && i < 100){
                prefix += "0" + std::to_string(i) + ".txt";
            }
            run_tpcase(prefix);
        }
    }
}