// overloading operators

template<typename T>
std::vector<T> operator-(const std::vector<T>& lhs,const std::vector<T>& rhs){
    if(lhs.size() != rhs.size()){std::cerr << "Wrong vector dimensions in operator - !" << std::endl; return std::vector<T>();}

    std::vector<T> result;
    for(int i=0;i < lhs.size();i++){result.push_back(lhs[i] - rhs[i]);}
    return result;
}
template<typename T>
std::vector<T> operator+(const std::vector<T>& lhs,const std::vector<T>& rhs){
    if(lhs.size() != rhs.size()){std::cerr << "Wrong vector dimensions in operator + !" << std::endl; return std::vector<T>();}

    std::vector<T> result;
    for(int i=0;i < lhs.size();i++){result.push_back(lhs[i] + rhs[i]);}
    return result;
}
template<typename T>
std::vector<T> operator+(const T lhs,const std::vector<T>& rhs){
    std::vector<T> result;
    for(int i=0;i < rhs.size();i++){result.push_back(rhs[i] + lhs);}
    return result;
}
template<typename T>
std::vector<T> operator+(const std::vector<T>& lhs,const T rhs){
    std::vector<T> result;
    for(int i=0;i < lhs.size();i++){result.push_back(lhs[i] + rhs);}
    return result;
}
template<typename T,typename G>
std::vector<T> operator/(const std::vector<T>& lhs,const G rhs){
    std::vector<T> result;
    for(int i=0;i < lhs.size();i++){result.push_back(lhs[i] / rhs);}
    return result;
}
template<typename T>
std::vector<T> operator/(const std::vector<T>& lhs,const std::vector<T>& rhs){
    if(lhs.size() != rhs.size()){std::cerr << "Wrong vector dimensions in operator / !" << std::endl; return std::vector<T>();}

    std::vector<T> result;
    for(int i=0;i < lhs.size();i++){result.push_back(lhs[i] / rhs[i]);}
    return result;
}
template<typename T,typename G>
std::vector<T> operator*(const std::vector<T>& lhs,const G rhs){
    std::vector<T> result;
    for(int i=0;i < lhs.size();i++){result.push_back(lhs[i] * rhs);}
    return result;
}
template<typename T,typename G>
std::vector<T> operator*(const G lhs,const std::vector<T>& rhs){
    std::vector<T> result;
    for(int i=0;i < rhs.size();i++){result.push_back(lhs * rhs[i]);}
    return result;
}
template<typename T>
std::vector<T> operator*(const std::vector<T>& lhs,const std::vector<T>& rhs){
    if(lhs.size() != rhs.size()){std::cerr << "Wrong vector dimensions in operator / !" << std::endl; return std::vector<T>();}

    std::vector<T> result;
    for(int i=0;i < lhs.size();i++){result.push_back(lhs[i] * rhs[i]);}
    return result;
}

// output

/// @brief prints a vector containing arbitrary type for which the operator << is defined. When width > 0, the width of the items will be width. Overload for a vector containing vector of arbitrary type.
template<typename T>
void print(std::vector<T> vec,int width){
    std::cout << "[";
    if(vec.size() == 0){std::cout << "]";}
    for(int i=0;i<vec.size();i++){
        std::cout << std::setw(width) << vec[i];
        if(i==vec.size()-1){std::cout << "]";}else{std::cout << " ";}
    }
    std::cout << std::endl;
}

/// @brief prints a vector containing arbitrary type for which the operator << is defined. When width > 0, the width of the items will be width. Overload for a vector containing vector of arbitrary type.
template<typename T>
void print(std::vector<std::vector<T> > vec, int width){
    for(int i=0;i<vec.size();i++){
        print(vec[i],width);
    }
}

/// @brief writes a vector containing arbitrary type for which the operator << is defined to the given fstream. When width > 0, the width of the items will be width. Overload for a vector containing vector of arbitrary type.
template<typename T>
void write_to_file(std::vector<T> vec, std::fstream& file,int width,std::string begin,std::string end){
    file << begin;
    if(vec.size() == 0){file << "]";}
    for(int i=0;i<vec.size();i++){
        file << std::setw(width) << vec[i];
        if(i==vec.size()-1){file << end;}else{file << ",";}
    }
    file << std::endl;
}

/// @brief writes a vector containing arbitrary type for which the operator << is defined to the given fstream. When width > 0, the width of the items will be width. Overload for a vector containing vector of arbitrary type.
template<typename T>
void write_to_file(std::vector<std::vector<T> > vec, std::fstream& file, int width,std::string begin,std::string end){
    for(int i=0;i<vec.size();i++){
        write_to_file(vec[i],file,width,begin,end);
    }
}

// miscellaneous

/// @brief The signum function. For x = 0, we have sgn(0) = 0. 
template<typename T>
int sgn(T x){return (T(0) < x) - (x < T(0));}

