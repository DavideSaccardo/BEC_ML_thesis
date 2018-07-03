#include "randomuniform.h"
#include <iostream>
#include <cassert>
#include "../Math/random.h"
#include "../system.h"
#include <cmath>
#include <random>
#include <vector>

using std::cout;
using std::endl;

using namespace std;

RandomUniform::RandomUniform(System* system, int numberOfParticles, int numberOfDimensions, int numberOfHiddenNodes, int numberOfVisibleNodes, double sigma,vector<double>& m_X, vector<double>& m_a_bias, vector<double>& m_b_bias, vector<std::vector<double>>& m_w,  double interactionSize, double timeStep, int numberOfParameters)   :
        InitialState(system) {
    assert(numberOfHiddenNodes > 0 && numberOfVisibleNodes > 0);

    m_system->setSigma                            (sigma);
    m_system->setSigma_squared                    (sigma*sigma);
    m_system->setNumberOfHiddenNodes              (numberOfHiddenNodes);
    m_system->setNumberOfVisibleNodes             (numberOfVisibleNodes);
    m_system->setNumberOfParticles                (numberOfParticles);
    m_system->setNumberOfDimensions               (numberOfDimensions);
    m_system->setinteractionSize                  (interactionSize);
    m_system->setTimeStep                         (timeStep);
    m_system->setSqrtTimeStep                     (sqrt(timeStep));


    setupInitialState(m_X, m_a_bias, m_b_bias, m_w);
}

void RandomUniform::setupInitialState(vector<double>& m_X, vector<double>& m_a_bias, vector<double>& m_b_bias, vector<std::vector<double>>& m_w) {

    //Set Sigma of the Gaussian distribution used to generate data
    double sigma_0=0.1; //0.5 for non interacting, 0.1 for interacting;

    int M = m_system->getNumberOfVisibleNodes();
    int N= m_system->getNumberOfHiddenNodes();



    for (int i=0; i < M; i++) {
        m_X[i]=(Random::nextDouble()-0.5);
        //cout<<"X: "<< m_X[i]<<endl;
    }

    m_system->computematrixdistance(m_X);


    for(int i=0; i<M; i++){
        m_a_bias[i]=Random::nextGaussian(0,sigma_0);
        //cout<<"a_bias: "<< m_a_bias[i]<<endl;
    }

    for(int i=0; i<N; i++){
        m_b_bias[i]=Random::nextGaussian(0,sigma_0);
        //cout<<"b_bias: "<< m_b_bias[i]<<endl;
    }

    for(int i=0; i<M; i++){
        for(int j=0;j<N; j++){
            m_w[i][j]=Random::nextGaussian(0,sigma_0);
        //cout<<"w: ["<<i<<" "<<j<<"]"<< m_w[i][j]<<endl;
        }
    }

}




