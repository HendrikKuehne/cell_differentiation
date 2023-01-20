/// @brief Structure representing one single cell
struct Cell{
    std::vector<double> p,m;

    std::vector<std::vector<double> > p_timeline,m_timeline;

    /// @brief Updates the protein and mRNA concentrations. The update is calculated using the last values in the timelines for m and p.
    /// @param Penv Protein concentration in the medium.
    /// @param J The matrix representation of the Gene Regulatory Network (GRN).
    void update(std::vector<double> Penv,std::vector<std::vector<double> > J){
        std::vector<double> old_p = *p_timeline.rbegin(), old_m = *m_timeline.rbegin();
        // calculating the change in protein concentration p and mRNA concentration m and adding it to p and m, respectively
        p = old_p + ((old_m - old_p) + D * (Penv - old_p)) * dt;
        m = old_m + gamma * (F(J,old_p) - old_m) * dt;

        // updating the timelines
        p_timeline.push_back(p);
        m_timeline.push_back(m);
    }

    Cell(){
        p = std::vector<double>();
        m = std::vector<double>();
        for(int i=0;i<theta.size();i++){p.push_back(cont_dist(engine)); m.push_back(cont_dist(engine));}

        // initializing the timelines
        p_timeline = std::vector<std::vector<double> >(); p_timeline.push_back(p);
        m_timeline = std::vector<std::vector<double> >(); m_timeline.push_back(m);
    }

    Cell(const Cell& oldcell,std::vector<double> noise){
        p = oldcell.p; m = oldcell.p;

        // cell reproduction is noisy
        p = p + noise;

        // initializing the timelines; the new cell inherits the timelines from the cell it emerged from
        p_timeline = oldcell.p_timeline; p_timeline[p_timeline.size()-1] = p;
        m_timeline = oldcell.m_timeline; m_timeline[m_timeline.size()-1] = m;
    }
};

/// @brief Matrix-vector multiplication
template<typename T>
std::vector<T> matvecmul(std::vector<std::vector<T> > matrix, std::vector<T> vector){
    int resdim = matrix.size();
    int indim = matrix[0].size();

    // checking dimensions
    for(auto itr = matrix.begin();itr < matrix.end();itr++){
        if(itr->size() != indim){
            std::cerr << "Wrong dimensions in matvecmul()!" << std::endl;
            return std::vector<T>();
        }
    }
    if(vector.size() != indim){
        std::cerr << "Wrong dimensions in matvecmul()!" << std::endl;
        return std::vector<T>();
    }

    // actual multiplication
    std::vector<T> result;
    for(int i=0;i<resdim;i++){
        result.push_back(0);
        for(int j=0;j<indim;j++){
            result[i] += matrix[i][j] * vector[j];
        }
    }

    return result;
}

/// @brief Sigmoidal activaion function f regulating the time evolution of mRNA
template<typename T>
T f(T x){return 1/(1 + std::exp(-beta*x));}
/// @brief Sigmoidal activaion function f regulating the time evolution of mRNA
std::vector<double> f(std::vector<double> x){
    for(std::vector<double>::iterator itr = x.begin();itr < x.end();itr++){*itr = 1/(1 + std::exp((-1) * beta * (*itr)));}
    return x;
}

/// @brief Gene expression level Function F for one cell.
/// @param J Representation of the GRN (gene regulatory network). If J[i][j] == 1, protein j activates gene i.
/// @param p expression levels of the proteins.
std::vector<double> F(std::vector<std::vector<double> > J, std::vector<double> p){return f(matvecmul(J,p) - theta);}

/// @brief Concentration of protein in the medium
/// @param cells A vector containing all the cells whose average protein concentration is being computed
std::vector<double> P(std::vector<Cell> cells){
    std::vector<double> result(theta.size(),0);
    for(std::vector<Cell>::iterator itr = cells.begin();itr < cells.end();itr++){result = result + itr->p;}
    result = result / cells.size();

    return result;
}

/// @brief Dividing all cells in the input
/// @param cells A vector containing all the old cells
std::vector<Cell> Mitosis(std::vector<Cell> cells){
    std::vector<Cell> result;
    std::vector<double> noise(5,0);

    // creating new cells from old ones
    for(std::vector<Cell>::iterator itr = cells.begin();itr < cells.end();itr++){
        // determining the noise - now woth a discrete random variable (like Julian did it) rather than a continuos one
        for(int i=0;i<theta.size();i++){
            // noise[i] = (2 * cont_dist(engine) - 1) * sigma;
            noise[i] = (2 * disc_dist(engine) - 1) * sigma;
        }

        // creating the new cells
        result.push_back(Cell(*itr,noise));
        result.push_back(Cell(*itr,(-1)*noise));
    }

    return result;
}

/// @brief Creates folder named after the current time and, within, it, creates a file parameters.txt in which the parameters of the simuation are stored
/// @param tdiv Time between cell divisions.
/// @param Nmax Maximum number of cells which we are simulating.
/// @param J The GRN we used for the simulation.
std::string make_documentation(double tdiv, int Nmax, std::vector<std::vector<double> > J){
    // exracting the time
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream dirname;
    dirname << "simulation_output/" << std::put_time(&tm, "%Y_%m_%d-%H_%M");
    std::ostringstream filename;
    filename << "simulation_output/" << std::put_time(&tm, "%Y_%m_%d-%H_%M") << "/abstract.txt";

    // cerating the directory and the file
    std::filesystem::create_directory(dirname.str());
    std::fstream file;
    file.open(filename.str(),std::ios::out);

    file << "Simulation parameters:\n    tdiv  = " << tdiv <<
                                  "\n    Nmax  = " << Nmax <<
                                  "\n    beta  = " << beta <<
                                  "\n    gamma = " << gamma <<
                                  "\n    dt    = " << dt <<
                                  "\n    D     = ";
    write_to_file(D,file);
    file << "    sigma = " << sigma << "\n    theta = ";
    write_to_file(theta,file);
    file << "\nGene Regulatory Network:\n";

    write_to_file(J,file);

    file.close();

    return dirname.str();
}

/// @brief Writes the timelines of all cells in cells to .txt files.
/// @param dirname Directory in which the files m_cell<index>.txt and p_cell<index>.txt are created for every index.
/// @param cells Vcetor containing all cells.
void write_timelines(std::string dirname, std::vector<Cell> cells){
    // iterating over all cells
    for(int i=0;i < cells.size();i++){
        // writing the mRNA timeline
        std::ostringstream filename;
        filename << dirname << "/m_cell" << i << ".txt";
        std::fstream file;
        file.open(filename.str(),std::ios::out);
        write_to_file(cells[i].m_timeline,file,0,"","");
        file.close();

        // writing the protein timeline
        filename = std::ostringstream();
        filename << dirname << "/p_cell" << i << ".txt";
        file.open(filename.str(),std::ios::out);
        write_to_file(cells[i].p_timeline,file,0,"","");
        file.close();
    }
}

/// @brief Checks if J is defined according to the paper, e.g. if it has the right dimensions, if all it's components are either 0,1 or -1 and if there are 10 non-zero components
/// @param J The GRN we used for simulation
bool check_J(std::vector<std::vector<double> > J){
    // checking the first dimension:
    if(not (J.size() == theta.size())){return false;}

    // checking the second dimension and the components:
    for(std::vector<std::vector<double> >::iterator itr1 = J.begin();itr1 < J.end();itr1++){
        if(itr1->size() != theta.size()){return false;}
        for(std::vector<double>::iterator itr2 = itr1->begin();itr2 < itr1->end();itr2++){
            if(abs(*itr2) != 1 and *itr2 != 0){return false;}
        }
    }

    return true;
}
