#include "system.h"
#include <cassert>
#include <cmath>
#include <vector>
#include "sampler.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"
#include <iostream>
#include <fstream>
#include <time.h>
#//include "conjugategradient.h"

using namespace std;

bool System::metropolisStepBruteForce(double GibbsValue, bool interaction,vector<double> &X, vector<double> a_bias, vector<double> b_bias, vector<std::vector<double>> w) { //Brute Force Metropolis method
    //Function to use Brute Force Metropolis
    int randparticle=Random::nextInt(getNumberOfVisibleNodes());

    vector <double> X_old (getNumberOfVisibleNodes());
    vector <double> X_new (getNumberOfVisibleNodes());

    X_old = X;
    X_new = X;

    vector <vector<double>> OldDistanceMatrix;
    OldDistanceMatrix=getDistanceMatrix();

    //choose a new move
    //int init =randparticle*m_numberOfDimensions;

    X_new[randparticle] = X_old[randparticle] + m_stepLength *( Random::nextDouble() - 0.5 );

    setDistanceMatrix (computematrixdistance(X_new));

    double psi_new = m_waveFunction->evaluate(GibbsValue, X_new, a_bias, b_bias, w);

    //   cout<<"new"<<psi_new<<endl;

    double prob = psi_new * psi_new / ( m_psiOld * m_psiOld );

    if (Random::nextDouble() < prob||1.0<prob){ // Accept the new move

        m_psiOld = psi_new;
        X        = X_new;

        return true;

    } else { // Don't accept accept the new move

        X = X_old;
        setDistanceMatrix(OldDistanceMatrix);

        return false;

    }
}


bool System::metropolisStepImportance(double GibbsValue, bool interaction,vector<double> &X, vector<double> a_bias, vector<double> b_bias, vector<std::vector<double>> w) { //Importance Sampling method
    //Function to use Importance Sampling

    int randparticle = Random::nextInt (getNumberOfVisibleNodes());

    vector <double> X_old (getNumberOfVisibleNodes());
    vector <double> X_new (getNumberOfVisibleNodes());
    X_old = X;
    X_new = X;

    vector <double> QF_old (m_numberOfVisibleNodes);
    vector <double> QF_new (m_numberOfVisibleNodes);
    QF_old = getQuantumForce();

    vector <vector<double>> OldDistanceMatrix;
    OldDistanceMatrix=getDistanceMatrix();

    X_new[randparticle] = X_old[randparticle] + 0.5 * m_QuantumForce[randparticle] * m_timeStep + m_sqrtTimeStep * ( Random::nextGaussian(0.0, 1.0) );

    setQuantumForce (m_waveFunction->QuantumForce(GibbsValue,X_new,a_bias,b_bias,w)); //UPDATE
    QF_new = m_QuantumForce;

    updateDistanceMatrix(X_new, randparticle);


    double GreensFunction = 0.0;

    GreensFunction = 0.5* (QF_old[randparticle] + QF_new[randparticle])
            * (0.5 * 0.5 * m_timeStep * (+QF_old[randparticle] - QF_new[randparticle])
               - X_new[randparticle] + X_old[randparticle] );

    GreensFunction = exp( GreensFunction );

    double psi_new = m_waveFunction->evaluate(GibbsValue,X_new,a_bias,b_bias, w);
    double prob    = GreensFunction * psi_new * psi_new / (m_psiOld * m_psiOld);

    // Accept  new move
    if ( (Random::nextDouble() < prob) || ( 1.0<prob )){

        m_psiOld = psi_new;
        X        = X_new;

        return true;
    }

    // Don't accept new move
    else{

        X = X_old;

        setQuantumForce   (QF_old);
        updateDistanceMatrix(X_old,randparticle);

        return false;
    }
}

bool System::GibbsSampling(vector<double> &X, vector<double> a_bias, vector<double> b_bias, vector<std::vector<double>> w){
    //Function to use Gibbs Sampling

    int N = getNumberOfHiddenNodes();
    int M = getNumberOfVisibleNodes();

    vector<double> Hidden(N);

    double sum;
    double sum2;
    double mean;
    double argument;

    for(int j = 0; j < N; j++){
        sum = 0;

        for(int i = 0; i < M; i++){
            sum += X[i] * w[i][j] / getSigma_squared();
            // cout<<"old X["<<i<<"]"<<X[i]<<endl;
        }

        argument = exp( - b_bias[j] - sum);
        Hidden[j]=1/(1+argument);
    }

    for( int i = 0; i < M; i++){

        sum2 = 0;

        for( int j = 0; j < N; j++){

            sum2 += Hidden[j] * w[i][j];

        }

        mean = a_bias[i] + sum2;
        X[i] = Random::nextGaussian(mean,getSigma());
        // cout<<"new X["<<i<<"]"<<X[i]<<endl;

    }

    return true;
}

void System::runMetropolisSteps(string method, int cycle, vector<double>&Gradient,int numberOfMetropolisSteps,bool interaction, vector<double> X, vector<double> a_bias, vector<double> b_bias, vector<std::vector<double>> w) {
    //Principal function of the whole code. Here the Monte Carlo method is

    double GibbsValue = 1.0;
    if(method=="Gibbs") GibbsValue=0.5;

    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;

    getSampler()->setStepNumber(0);
    getSampler()->setAcceptedNumber(0);
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
    m_sampler->setDimensionOfGradient(m_numberOfParameters);

    // Initial values
    setDistanceMatrix(computematrixdistance(X));

    m_psiOld = m_waveFunction->evaluate(GibbsValue, X, a_bias, b_bias,w);
    setQuantumForce(m_waveFunction->QuantumForce(GibbsValue, X,a_bias,b_bias,w));

    vector<double> temp (m_numberOfParameters);
    setCumulativeEnGradient(temp);
    setCumulativeGradient(temp);
    setGradient(temp);
    setEnGradient_average(temp);

    bool acceptedStep;

    for (int i=0; i < numberOfMetropolisSteps; i++) {

        if( method == "MetropolisBruteForce" ) acceptedStep =  metropolisStepBruteForce(GibbsValue,interaction, X, a_bias, b_bias,w);   //run the system with Brute Force Metropolis
        if( method == "MetropolisImportance" ) acceptedStep =  metropolisStepImportance(GibbsValue, interaction,X, a_bias, b_bias,w);    //run the system with Importance Sampling
        if( method == "Gibbs"                ) acceptedStep =  GibbsSampling(X,a_bias,b_bias,w);

        m_sampler->sample(GibbsValue, acceptedStep,interaction,X,a_bias,b_bias,w);  //sample results and write them to file
        m_sampler->writeToFile(i);

    }
    m_sampler->closeDataFile();

    m_sampler->computeAverages(Gradient);
}


void System::StochasticGradientDescent(vector<double> Gradient, vector<double> X, vector<double>& a_bias, vector<double>& b_bias, vector<std::vector<double>>& w){
    //SGD: parameters are updated according to the Gradient of energy
    for( int i = 0; i < m_numberOfVisibleNodes; i++){

        a_bias[i] -= m_learningRate * Gradient[i];
        // cout<<"a:["<<i<<"]"<<a_bias[i]<<endl;

    }

    for(int i = 0; i < m_numberOfHiddenNodes; i++){

        b_bias[i] -= m_learningRate * Gradient[i + m_numberOfVisibleNodes];
        //cout<<"b:["<<i<<"]"<<b_bias[i]<<endl;

    }

    int z = m_numberOfVisibleNodes + m_numberOfHiddenNodes;

    for(int i = 0; i < m_numberOfVisibleNodes; i++){

        for(int j = 0; j < m_numberOfHiddenNodes; j++){

            w[i][j] -= m_learningRate * Gradient[z];

            z++;
            //      cout<<"w:["<<i<<"]["<<j<<"]"<< w[i][j]<<endl;

        }
    }
}

int System::getNumberOfParameters() const
{
    return m_numberOfParameters;
}

void System::setNumberOfParameters(int numberOfParameters)
{
    m_numberOfParameters = numberOfParameters;
}


vector<double> System::GradientParameters(double GibbsValue, vector<double> X, vector<double>& a_bias, vector<double>& b_bias, vector<std::vector<double>>& w){
    //Return Gradient of wf with respect to the parameters

    vector<double> GradientPsi (m_numberOfParameters);
    vector<double> argument    (m_numberOfHiddenNodes);

    double sum;

    for(int j = 0; j < m_numberOfHiddenNodes; j++){

        sum=0;

        for(int i = 0; i < m_numberOfVisibleNodes; i++){

            sum += X[i] * w[i][j] / getSigma_squared();

        }

        argument[j] = b_bias[j] + sum;
    }

    for(int i = 0; i < m_numberOfVisibleNodes; i++){
        GradientPsi[i] = ( X[i] - a_bias[i] ) * GibbsValue / getSigma_squared();
    }


    for(int i = m_numberOfVisibleNodes; i < m_numberOfHiddenNodes + m_numberOfVisibleNodes; i++){
        GradientPsi[i] = GibbsValue / ( 1 + exp ( - argument[i - m_numberOfVisibleNodes] ) );
    }

    int z = m_numberOfVisibleNodes + m_numberOfHiddenNodes;
    for(int i = 0; i < m_numberOfVisibleNodes; i++){

        for(int j=0; j<m_numberOfHiddenNodes; j++){

            GradientPsi[z] = GibbsValue * X[i] / ( getSigma_squared() * ( 1 + exp ( - argument[j] ) ) );
            z++;

        }

    }

    return GradientPsi;
}

std::vector<double> System::getCumulativeGradient() const
{
    return m_cumulativeGradient;
}

void System::setCumulativeGradient(const std::vector<double> &cumulativeGradient)
{
    m_cumulativeGradient = cumulativeGradient;
}

std::vector<double> System::getCumulativeEnGradient() const
{
    return m_cumulativeEnGradient;
}

void System::setCumulativeEnGradient(const std::vector<double> &cumulativeEnGradient)
{
    m_cumulativeEnGradient = cumulativeEnGradient;
}

std::vector<double> System::getGradient() const
{
    return m_Gradient;
}

void System::setGradient(const std::vector<double> &Gradient)
{
    m_Gradient = Gradient;
}

std::vector<double> System::getEnGradient_average() const
{
    return m_EnGradient_average;
}

void System::setEnGradient_average(const std::vector<double> &EnGradient_average)
{
    m_EnGradient_average = EnGradient_average;
}


int System::getNumberOfVisibleNodes() const
{
    return m_numberOfVisibleNodes;
}

void System::setNumberOfVisibleNodes(int numberOfVisibleNodes)
{
    m_numberOfVisibleNodes = numberOfVisibleNodes;
}

int System::getNumberOfHiddenNodes() const
{
    return m_numberOfHiddenNodes;
}

void System::setNumberOfHiddenNodes(int numberOfHiddenNodes)
{
    m_numberOfHiddenNodes = numberOfHiddenNodes;
}

double System::getSigma() const
{
    return m_sigma;
}

void System::setSigma(double sigma)
{
    m_sigma = sigma;
}

double System::getSigma_squared() const
{
    return m_sigma_squared;
}

void System::setSigma_squared(double sigma_squared)
{
    m_sigma_squared = sigma_squared;
}

double System::getLearningRate() const
{
    return m_learningRate;
}

void System::setLearningRate(double learningRate)
{
    m_learningRate = learningRate;
}


void System::printOut(int cycle, int totalNumberOfCycles)
{
    // m_sampler->computeAverages();

    m_sampler->printOutputToTerminal(cycle);
    //m_sampler->writeToFile(cycle);

    //    if(cycle==totalNumberOfCycles) {
    //        m_sampler->closeDataFile();
    //    }

}



//double System::computedistance(int i){
//    double temp=0;
//    for(int j=0;j<m_numberOfDimensions;j++){
//        temp+=m_particles.at(i).getPosition()[j]*m_particles.at(i).getPosition()[j];
//    }
//    return sqrt(temp);
//}

int System::computeIndex(int index){
    int init=index;
    //cout<<"begin"<<index<<endl;
    if((index%m_numberOfDimensions)==0) {
        //cout<<"end"<<init<<endl;
        return index;
    }

    while(index>m_numberOfDimensions){
        index -=m_numberOfDimensions;
    }
    init -=index;
    // cout<<"end"<<init<<endl;
    return init;
}

void System::updateDistanceMatrix(std::vector<double> m_X, int randparticle){

    double temp = 0;
    double  init=computeIndex(randparticle);
    //cout<<"ehi"<<endl;
    double part=init/m_numberOfDimensions;

    //int init=randparticle*m_numberOfDimensions;
    for (int j = 0; j < init; j+=m_numberOfDimensions){

        temp = 0;

        for (int d = 0; d < m_numberOfDimensions; d++){

            temp += (m_X[j+d] - m_X[init+d]) *
                    (m_X[j+d] - m_X[init+d]);
            // cout<<temp<<endl;

        }
        m_distanceMatrix[part][j/m_numberOfDimensions] = sqrt(temp);

        // cout<<part<<endl;
        //cout<<j/m_numberOfDimensions<<endl;
        //        if (m_distanceMatrix[randparticle][j] < getinteractionSize()){
        //            return true;
        //        }

        m_distanceMatrix[j/m_numberOfDimensions][part] = m_distanceMatrix[part][j/m_numberOfDimensions];
    }

    for (int j = init+m_numberOfDimensions; j < m_numberOfVisibleNodes; j+=m_numberOfDimensions){

        temp = 0;

        //   cout<<"ehi"<<endl;

        for (int d = 0; d < m_numberOfDimensions; d++){

            temp += (m_X[j+d] - m_X[init+d]) *
                    (m_X[j+d] - m_X[init+d]);

        }

        m_distanceMatrix[part][j/m_numberOfDimensions] = sqrt(temp);
        //        if (m_distanceMatrix[randparticle][j] < getinteractionSize()){
        //            return true;
        //        }
        m_distanceMatrix[j/m_numberOfDimensions][part] = m_distanceMatrix[part][j/m_numberOfDimensions];

    }

}

std::vector<vector<double>> System::computematrixdistance(vector<double>& m_X){

    vector<vector<double>> distancematrix(m_numberOfParticles, vector<double>(m_numberOfParticles));
    double temp=0;
    int j=0;
    int z;
    int k=0;

    while(j < m_numberOfVisibleNodes){

        temp = 0;
        z    = 0;

        for(int i = 0;i < j; i += m_numberOfDimensions){

            for(int q = 0; q < m_numberOfDimensions; q++){

                temp += ( m_X[i + q] - m_X[j + q] ) * ( m_X[i + q] - m_X[j + q] );

            }
            //temp=(m_X[i]-m_X[j])*(m_X[i]-m_X[j])+(m_X[i+1]-m_X[j+1])*(m_X[i+1]-m_X[j+1])+(m_X[i+2]-m_X[j+2])*(m_X[i+2]-m_X[j+2]);

            distancematrix[z][k] = sqrt(temp);
            //  cout<<distancematrix[z][k]<<endl;
            distancematrix[k][z] = distancematrix[z][k];

            z++;
        }

        j += m_numberOfDimensions;
        k++;
    }
    return distancematrix;
}

//double System::computedistanceABS(int i, int j){
//    double temp=0;
//    for(int k=0;k<m_numberOfDimensions;k++){
//        temp+=(m_particles.at(i).getPosition()[k] - m_particles.at(j).getPosition()[k]) *
//                (m_particles.at(i).getPosition()[k] - m_particles.at(j).getPosition()[k]);
//    }
//    return sqrt(temp);{
//    }
//}

void System::openDataFile(string filename){
    m_sampler->openDataFile(filename);
}

double System::getinteractionSize() const
{
    return m_interactionSize;
}

void System::setinteractionSize(double interactionSize)
{
    m_interactionSize = interactionSize;
}

double System::getTimeStep() const
{
    return m_timeStep;
}

void System::setTimeStep(double timeStep)
{
    m_timeStep = timeStep;
}

double System::getSqrtTimeStep() const
{
    return m_sqrtTimeStep;
}

void System::setSqrtTimeStep(double sqrtTimeStep)
{
    m_sqrtTimeStep = sqrtTimeStep;
}

std::vector<vector<double> > System::getDistanceMatrix() const
{
    return m_distanceMatrix;
}

double System::getDistanceMatrixij(int i, int j) const
{
    return m_distanceMatrix[i][j];
}

double System::getPsiOld() const
{
    return m_psiOld;
}

void System::setPsiOld(double psiOld)
{
    m_psiOld = psiOld;
}

std::vector<double> System::getQuantumForce() const
{
    return m_QuantumForce;
}

void System::setQuantumForce(const std::vector<double> &QuantumForce)
{
    m_QuantumForce = QuantumForce;
}

void System::setDistanceMatrix(const std::vector<vector<double> > &distanceMatrix)
{
    m_distanceMatrix = distanceMatrix;
}


void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
}


double System::findEnergyDerivative()
{

    double meanEnergy      =  getSampler()->getCumulativeEnergy()          / (m_numberOfMetropolisSteps*getEquilibrationFraction());
    double meanWFderiv     =  getSampler()->getCumulativeWFderiv()         / (m_numberOfMetropolisSteps*getEquilibrationFraction());
    double meanWFderivEloc =  getSampler()->getCumulativeWFderivMultEloc() / (m_numberOfMetropolisSteps*getEquilibrationFraction());

    return 2 * (meanWFderivEloc - meanEnergy*meanWFderiv);
}

