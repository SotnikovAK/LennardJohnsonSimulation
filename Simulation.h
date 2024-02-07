#pragma once

#include "Element.h"

class Simulation
{
private:



    const int WindowSize[2] = { 1000, 1000 };
    const int SpawnBorderEpsillon[2] = { 100, 100 };

    const float dt = 0.001;
    int n;

    float BorderSize;
    float landscape;
    int flag_animation;

    double U = 0.0;


    std::list<ElementarElement> ElementList;
    sf::RenderWindow window;

public:
    Simulation(float BorderSize, int n, int flag_animation);

    void WindowInteraction();

    void CreationElements(float BorderSize, int n);
};
