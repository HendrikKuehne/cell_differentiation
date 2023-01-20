#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<sstream>
#include<filesystem>
#include<random>
#include<chrono>

// global parameters
extern const double beta,gamma,sigma;
extern double dt;
extern const std::vector<double> theta,D;

// random number generation
extern std::default_random_engine engine;
extern std::uniform_real_distribution<double> cont_dist;
extern std::discrete_distribution<double> disc_dist;

// operator overloads
template<typename T>
std::vector<T> operator-(const std::vector<T>& lhs,const std::vector<T>& rhs);
template<typename T>
std::vector<T> operator+(const std::vector<T>& lhs,const std::vector<T>& rhs);
template<typename T>
std::vector<T> operator+(const T lhs,const std::vector<T>& rhs);
template<typename T>
std::vector<T> operator+(const std::vector<T>& lhs,const T rhs);
template<typename T,typename G>
std::vector<T> operator/(const std::vector<T>& lhs,const G rhs);
template<typename T>
std::vector<T> operator/(const std::vector<T>& lhs,const std::vector<T>& rhs);
template<typename T,typename G>
std::vector<T> operator*(const std::vector<T>& lhs,const G rhs);
template<typename T,typename G>
std::vector<T> operator*(const G lhs,const std::vector<T>& rhs);
template<typename T>
std::vector<T> operator*(const std::vector<T>& lhs,const std::vector<T>& rhs);

// structs

/// @brief Structure representing one single cell
struct Cell;

// functions
template<typename T>
std::vector<T> matvecmul(std::vector<std::vector<T> > matrix, std::vector<T> vector);
template<typename T>
T f(T x);
std::vector<double> f(std::vector<double> x);
std::vector<double> F(std::vector<std::vector<double> > J, std::vector<double> p);
std::vector<double> P(std::vector<Cell> cells);
std::vector<Cell> Mitosis(std::vector<Cell> cells);
std::string make_documentation(double tdiv, int Nmax, std::vector<std::vector<double> > J);
void write_timelines(std::string dirname, std::vector<Cell> cells);
bool check_J(std::vector<std::vector<double> > J);

// output

template<typename T>
void print(std::vector<T> vec,int width = 0);
template<typename T>
void print(std::vector<std::vector<T> > vec, int width = 0);
template<typename T>
void write_to_file(std::vector<T> vec, std::fstream& file,int width = 0,std::string begin = "[",std::string end = "]");
template<typename T>
void write_to_file(std::vector<std::vector<T> > vec, std::fstream& file, int width = 0,std::string begin = "[",std::string end = "]");

// miscellaneous

template<typename T>
int sgn(T x);