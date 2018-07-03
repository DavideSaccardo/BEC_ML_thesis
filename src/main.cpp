#include <iostream>
#include <random>
#include <cmath>
#include "system.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"
#include <chrono>
#include <string>
#include <vector>
#include <fstream>

using namespace std;
//std::ofstream ofile;


int main(){

    int numberOfParticles   = 1;        //Number of particles of the system considered
    int numberOfDimensions  = 3;        //Number of dimensions
    int numberOfHiddenNodes = 3;        //Number of hidden nodes
    int numberOfVisibleNodes =numberOfDimensions*numberOfParticles;

    double sigma =1.0;

    int numberOfParameters = numberOfVisibleNodes + numberOfHiddenNodes
            + numberOfVisibleNodes * numberOfHiddenNodes;


    bool interaction = true; //true = interaction, false = no interaction


    std::vector<double> X(numberOfVisibleNodes);

    std::vector<double> a_bias=std::vector<double>(numberOfVisibleNodes);

    std::vector<double> b_bias=std::vector<double>(numberOfHiddenNodes);

    std::vector<std::vector<double>> w(numberOfVisibleNodes, vector<double>(numberOfHiddenNodes));


    double beta                = 2.82843;                               // for interacting case: beta=2.82843
    int    numberOfSteps       = (int) 100000;                          // Number of Monte Carlo steps: 300000 for non-interacting, 100000 for interacting
    double interactionSize     = 0.0043;                                // for interacting case: interactionSize=0.0043;

    double timeStep            = 0.5;                                   // Importance sampling time step
    double stepLength          = 1.0;                                   // Metropolis step length.
    double omega               = 1.0;                                   // Oscillator frequency.
    double omega_z             = beta;                                  // Oscillator frequency in z-direction
    double equilibration       = 0.2;                                   // Amount of the total steps used for equilibration.

    double learning_rate       = 0.03;
    int    TotalNumberOfCycles = 1200;                                  // Total number of SGD cycles

    //Uncomment the suitable string to choose the desired sampling method
    //string method = "MetropolisBruteForce";
    string method = "MetropolisImportance";
    // string method = "Gibbs";

    // string filename = "0";
    // string filename = "Learn" + to_string(learning_rate) + "_TimeStep_"+to_string(timeStep) + "_Hidden_" + to_string(numberOfHiddenNodes)  + ".dat";
    // Set filename to "0" to stop from writing to file

    string filename_cycle_data="0";
    string finalFilename ="0";

    if( method == "MetropolisBruteForce" ) filename_cycle_data =  "bruCycleDataI_" + to_string(stepLength) + "_n_" + to_string(learning_rate) +  "_Np_" + to_string(numberOfParticles) + "_Nd_" + to_string(numberOfDimensions) + "_NH_" + to_string(numberOfHiddenNodes) +  "_w_" + to_string(omega) +  ".dat";
    if( method == "MetropolisImportance" ) filename_cycle_data = "impCycleDataI_" + to_string(timeStep) + "_n_" + to_string(learning_rate)+ "_Np_" + to_string(numberOfParticles) + "_Nd_" + to_string(numberOfDimensions) + "_NH_" + to_string(numberOfHiddenNodes) +  "_w_" + to_string(omega) +  ".dat";
    if( method == "Gibbs"                ) filename_cycle_data =  "gibCycledataI2_s_" + to_string(sigma) +"_n_" + to_string(learning_rate)+  "_Np_" + to_string(numberOfParticles) + "_Nd_" + to_string(numberOfDimensions) + "_NH_" + to_string(numberOfHiddenNodes) +  "_w_" + to_string(omega) +  ".dat";

    // Instantaneous energy data file (of bigger run after SGD)
    if( method == "MetropolisBruteForce" ) finalFilename =  "finalBruCycleDataI_" + to_string(stepLength) + "_n_" + to_string(learning_rate)+ "_Np_" + to_string(numberOfParticles) + "_Nd_" + to_string(numberOfDimensions) + "_NH_" + to_string(numberOfHiddenNodes) +  "_w_" + to_string(omega) +  ".dat";
    if( method == "MetropolisImportance" ) finalFilename = "finalImpCycleDataI_" + to_string(timeStep) + "_n_" + to_string(learning_rate)+ "_Np_" + to_string(numberOfParticles) + "_Nd_" + to_string(numberOfDimensions) + "_NH_" + to_string(numberOfHiddenNodes) +  "_w_" + to_string(omega) +  ".dat";
    if( method == "Gibbs"                ) finalFilename =  "finalGibCycledataI_s_" + to_string(sigma) +"_n_" + to_string(learning_rate)+ "_Np_" + to_string(numberOfParticles) + "_Nd_" + to_string(numberOfDimensions) + "_NH_" + to_string(numberOfHiddenNodes) +  "_w_" + to_string(omega) +  ".dat";

    //string filename = "parameter_Np_"+to_string(numberOfParticles)+"_Nd_"+to_string(numberOfDimensions)+"_Hidden_"+to_string(numberOfHiddenNodes);
    //ofile.open(filename);

    System* system = new System();
    system->setHamiltonian              (new HarmonicOscillator(system, omega, omega_z));
    system->setInitialState             (new RandomUniform(system, numberOfParticles, numberOfDimensions, numberOfHiddenNodes, numberOfVisibleNodes, sigma, X, a_bias, b_bias, w, interactionSize, timeStep, numberOfParameters));
    system->setWaveFunction             (new SimpleGaussian(system));
    //system->openDataFile              (filename_cycle_data);
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);



    vector <double> Gradient(numberOfParameters);

    //Compute computational time

    auto start = std::chrono::system_clock::now();

    string filename_cycle_energy="0";
    for(int cycles = 0; cycles < TotalNumberOfCycles; cycles ++){

        filename_cycle_energy="Np_"+to_string(numberOfParticles)+"_Nd_"+to_string(numberOfDimensions)+"_Hidden_"+to_string(numberOfHiddenNodes)+"_cycle_"+to_string(cycles)+".dat"; //"_Importance_learn_"+to_string(learning_rate)+"_cycle_"+to_string(cycles)+".dat";

        system->openDataFile              (filename_cycle_energy);

        system->setLearningRate           (learning_rate);

        system->setNumberOfParameters     (numberOfParameters);

        system->runMetropolisSteps        (method, cycles,Gradient,numberOfSteps,interaction, X, a_bias, b_bias, w);

        system->StochasticGradientDescent (Gradient,X,a_bias,b_bias,w);

        system->printOut                  (cycles,TotalNumberOfCycles);

    }

    if (interaction == false){

        int finalNumberOfSteps = 10000000;

        filename_cycle_energy="Np_"+to_string(numberOfParticles)+"_Nd_"+to_string(numberOfDimensions)+"_Hidden_"+to_string(numberOfHiddenNodes)+"_Importance.dat"; //_timestep_"+to_string(learning_rate)+"_cycle_"+to_string(TotalNumberOfCycles)+".dat";

        system->openDataFile              (filename_cycle_energy);
        system->runMetropolisSteps        (method, TotalNumberOfCycles, Gradient,finalNumberOfSteps,interaction, X, a_bias, b_bias, w);
        system->printOut                  (TotalNumberOfCycles, TotalNumberOfCycles);
        //system->writeToFile             (X,a_bias,b_bias,w);
    }

    //time

    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> diff = end - start;

    std::cout << " Computation time = " << diff.count() / 60.0 << " min\n" << endl; //display run time

    return 0;
}


