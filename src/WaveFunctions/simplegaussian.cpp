#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include <vector>
#include <algorithm>
#include "wavefunction.h"
#include "../system.h"
#include <iostream>

using namespace std;
using std::vector;

//parameters used to fix ML in the elliptic harmonic potential case
double beta = pow(2.82843,1.0/4.0);
double beta2 = 2.82843*2.82843;


SimpleGaussian::SimpleGaussian(System* system)://, vector<double> a_bias, vector<double> b_bias, vector<vector<double>> w) :
    WaveFunction(system) {
    m_numberOfParameters = 3;
    m_parameters.reserve(3);

}

double SimpleGaussian::evaluate(double GibbsValue,vector<double>& X, vector<double>& a_bias, vector<double>& b_bias, vector<std::vector<double>>& w) {

    double first_sum = 0;
    double prod = 1;
    int N = m_system->getNumberOfHiddenNodes();
    int M = m_system->getNumberOfVisibleNodes();

    for(int j = 0; j < m_system->getNumberOfParticles(); j++){

        for(int i = 0; i < j; i++){

            //cout<<m_system->getDistanceMatrixij(i,j)<<endl;

            if(m_system->getDistanceMatrixij(i,j)<=m_system->getinteractionSize()) {

                cout<<"Too close"<<endl;

                return 0;}

        }
    }


    for (int i = 0; i < M; i++){

        if((i+1)%3==0) {

            first_sum += (X[i]*beta-a_bias[i])*(X[i]*beta-a_bias[i]);

        }   else first_sum += (X[i]-a_bias[i])*(X[i]-a_bias[i]);
    }

    first_sum /= 2*m_system->getSigma_squared();

    first_sum = exp(-first_sum*GibbsValue);

    for (int j = 0; j < N; j++){

        double second_sum = 0;

        for (int i = 0; i < M; i++){

            if((i+1)%3==0) second_sum += beta*X[i]*w[i][j]*beta2;

            else second_sum += X[i]*w[i][j];
        }

        second_sum /= m_system->getSigma_squared();

        prod *= 1 + exp(b_bias[j] + second_sum);
    }

    if(GibbsValue == 0.5) { return first_sum*sqrt(prod);}

    return first_sum*prod;
}


double SimpleGaussian::computeDoubleDerivative(double GibbsValue, vector<double>& X, vector<double>& a_bias, vector<double>& b_bias, vector<std::vector<double>>& w) {

    int N = m_system->getNumberOfHiddenNodes();
    int M = m_system->getNumberOfVisibleNodes();

    vector <double> argument(N);

    double firstsum  = 0.0;
    double secondsum = 0.0;
    double kinetic   = 0.0;

    double temp2;
    double temp3;
    double sum;

    for(int j = 0; j < N; j++){
        sum = 0;

        for(int i = 0; i < M; i++){
            if((i+1)%3==0)sum += X[i] * w[i][j] *beta*beta2/ m_system->getSigma_squared();
            else sum += X[i] * w[i][j] / m_system->getSigma_squared();
        }

        argument[j] = exp( - b_bias[j] - sum);
        // cout<<argument[j]<<endl;
    }

    for(int i = 0 ; i < M;i++){

        temp2 = 0;
        temp3 = 0;

        for(int j = 0; j < N; j++){

            //double temp4 = exp(- argument[j]);
            double expon = 1.0 + argument[j];

            if((i+1)%3==0) {

                temp2 += w[i][j] * 1.0*beta2 / (expon);

                temp3 += w[i][j] * w[i][j]*beta2*beta2 * argument[j] / (expon*expon);

            } else {

                temp2 += w[i][j] * 1.0 / (expon);

                temp3 += w[i][j] * w[i][j] * argument[j] / (expon*expon);

            }
        }

        if((i+1)%3==0) firstsum = - ( X[i]*beta - a_bias[i] ) + temp2;

        else firstsum = - ( X[i] - a_bias[i] ) + temp2;

        firstsum /= m_system->getSigma_squared();

        secondsum = temp3 / ( m_system->getSigma_squared() * m_system->getSigma_squared() )
                    - 1.0 / m_system->getSigma_squared();

        //cout<<secondsum<<endl;
        //cout<<"somma"<<secondsum<<endl;

        kinetic += firstsum * firstsum*GibbsValue*GibbsValue + secondsum*GibbsValue;
    }

    //cout<<"firstsum: "<<firstsum<<endl;
    //cout<<"secondsum: " <<secondsum<<endl;
    //cout<<"third: "<<-M/m_system->getSigma_squared()<<endl;
    return kinetic;
}


std::vector<double> SimpleGaussian::QuantumForce(double GibbsValue, vector<double>& X, vector<double>& a_bias, vector<double>& b_bias, vector<std::vector<double>>& w) {
    //Function to comput the Quantum Force for the Importance Sampling method

    //double temp2;
    int N = m_system->getNumberOfHiddenNodes();
    int M = m_system->getNumberOfVisibleNodes();

    vector <double> QuantumForce(M);
    vector <double> argument(N);
    vector <double> temp2(M);

    double sum;

    for (int j = 0; j < N; j++){

        sum = 0;

        for (int i = 0; i < M; i++){
            if((i+1)%3==0) sum +=beta* X[i] * w[i][j]*beta2 / m_system->getSigma_squared();
            else sum += X[i] * w[i][j] / m_system->getSigma_squared();
        }
        argument[j] = b_bias[j] + sum;
    }

    for(int i = 0; i < M; i++){

        temp2[i] = 0;

        for(int j = 0; j < N; j++){

            double temp4 = exp( - argument[j] );
            double expon = 1 + temp4;
            if((i+1)%3==0){
            temp2[i] += w[i][j]*beta2 / ( expon );
            } else {
                temp2[i] += w[i][j] / ( expon );
            }
        }

        if((i+1)%3==0) QuantumForce[i] = 2 * ( - (X[i]*beta - a_bias[i] ) + temp2[i] ) * GibbsValue / m_system->getSigma_squared();
        else QuantumForce[i] = 2 * ( - (X[i] - a_bias[i] ) + temp2[i] ) * GibbsValue / m_system->getSigma_squared();
        //cout<<QuantumForce[i]<<endl;
    }

    return QuantumForce;
}


//    double a = m_system->getinteractionSize() ;
//    double constant;
//    double R_kj;
//    double dimension=m_system->getNumberOfDimensions();
//    double number =m_system->getNumberOfParticles();
//    std::vector<std::vector<double>> QuantumForce(dimension,vector<double>(number));
//    for (int d = 0; d < m_system->getNumberOfDimensions(); d++){
//        for (int k = 0; k < m_system->getNumberOfParticles(); k++){
//    QuantumForce[d][k] = -2 * (m_parameters[d]*particles.at(k).getPosition()[d]);
//        for (int j = 0; j < k; j++){
//                R_kj = m_system->getDistanceMatrixij(k,j);
//                constant = 2*a / (R_kj*R_kj*(R_kj-a));

//                    QuantumForce[d][k] += (particles.at(k).getPosition()[d] - particles.at(j).getPosition()[d]) * constant;

//            }
//        }
//    }
//    return QuantumForce;
//}

