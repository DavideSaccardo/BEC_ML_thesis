#pragma once
#include <vector>
#include "../particle.h"

class InitialState {

public:

    InitialState(class System* system);
    virtual void setupInitialState(std::vector<double> &m_X, std::vector<double> &m_a_bias, std::vector<double> &m_b_bias, std::vector<std::vector<double> > &m_w) = 0;

protected:

    class System* m_system = nullptr;
    class Random* m_random = nullptr;

};

