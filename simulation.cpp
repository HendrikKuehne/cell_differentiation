#include"header.h"
#include"functions.h"
#include"infrastructure.h"

/*
 * The expression levels of the k genes and proteins will be stored in vectors (std::vector<double>)
 * such that Jij becomes a matrix (std::vector<std::vector<double> >). There are two vectors m and p
 * for each cell which contain the gene and protein expression levels, respectively.
 */

// global parameters
const double beta = 40;                                             // slope of the sigmoidal function f
const double gamma = 6;                                             // "rate of change in mRNA"
const std::vector<double> theta = {-0.01,-0.03,0.02,0.01,-0.02};    // gene expression thresholds; theta.size() determines the number of proteins and genes in the rest of the program
const std::vector<double> D = {0,0,0,0,0.4};                        // Diffusion constants
const double sigma = 1e-3;                                          // Noise level during cell division
double dt = .01;                                                    // infinitesimal time difference

// random number generator
std::default_random_engine engine;
std::uniform_real_distribution<double> cont_dist(0,1);
std::discrete_distribution<double> disc_dist({1,1});

int main(int argc, char *argv[]){
    // parameters of our simulation
    double tdiv = 500;
    int Nmax = 128;
    if(argc == 2){
        tdiv = std::stod(argv[1]);
    }else if(argc == 3){
        tdiv = std::stod(argv[1]);
        Nmax = std::stoi(argv[2]);
    }else if(argc == 4){
        tdiv = std::stod(argv[1]);
        Nmax = std::stoi(argv[2]);
        dt = std::stod(argv[3]);
    }else if(argc > 4){
        std::cerr << "Command-line arguments: ./simulation.app <tdiv> <Nmax> <dt>" << std::endl;
        return -1;
    }
    std::cout << "Simulating with\n  tdiv = " << tdiv << "\n  Nmax = " << Nmax << "\n  dt = " << dt << std::endl << std::endl;

    // defining the GRN that we'll use
    // std::vector<std::vector<double> > J{{ 0, 0, 0, 0, 0},   // Figure 1a
    //                                     {-1, 0, 0, 0, 0},
    //                                     { 0, 1, 0, 0, 0},
    //                                     { 1, 0,-1, 1,-1},
    //                                     { 1, 0,-1, 1,-1}};
    std::vector<std::vector<double> > J{{ 0, 0, 0, 0, 0},   // Figure 1b
                                        { 1, 0, 0, 0, 0},
                                        { 0, 1,-1, 0, 0},
                                        { 0,-1, 1, 0, 1},
                                        { 0, 1,-1,-1, 1}};
    // std::vector<std::vector<double> > J{{ 0, 0, 0, 0, 0},   // Figure 1c
    //                                     { 0, 0,-1, 0, 0},
    //                                     { 0, 0, 0,-1, 1},
    //                                     { 1, 0,-1, 0,-1},
    //                                     { 0, 1,-1, 1,-1}};
    // std::vector<std::vector<double> > J{{ 0, 1, 0, 0, 0},   // Figure 1d
    //                                     {-1, 0, 1, 0, 0},
    //                                     { 0,-1, 0, 1, 0},
    //                                     { 0, 0,-1, 1, 1},
    //                                     { 0,-1, 0,-1, 0}};

    if(not check_J(J)){
        std::cerr << "J is wrongly defined!" << std::endl;
        return -1;
    }

    // seed for the random number generator
    engine.seed(std::chrono::system_clock::now().time_since_epoch().count());

    // defining a vector which will contain all the cells
    std::vector<Cell> cells(1,Cell());

    // infrastructure
    double t = 0;
    std::vector<double> Penv;

    // carrying out the simulation
    while(cells.size() <= Nmax and t < tdiv){
        // updating mRNA and protein concentrations
        Penv = P(cells);

        for(auto itr = cells.begin();itr < cells.end();itr++){
            itr->update(Penv,J);
        }

        // time passed
        t++;

        // cell division, if enough time passed
        if(t >= tdiv and cells.size() < Nmax){
            std::cout << "Cell division from N = " << cells.size() << " to ";
            cells = Mitosis(cells);
            std::cout << "N = " << cells.size() << std::endl;
            t = 0;
        }
    }

    // creating the directory to store our data in
    std::string dirname = make_documentation(tdiv,Nmax,J);

    // writing the timelines of all the cells to files
    write_timelines(dirname,cells);

    return 0;
}