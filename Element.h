#pragma once
#pragma once

#include "libraries.h"
#include "Constants.h"
#include "Quaternion.h"

class ElementarElement
{
public:
    float M = MassParticale;

    float Vx = 0.0f;
    float Vy = 0.0f;
    float Vz = 0.0f;

    float Wx = 0.0f;
    float Wy = 0.0f;
    float Wz = 0.0f;

    float R = 0.1;

    float x = 0.0f;
    float y = 0.0f;
    float z = 0.0f;

    float x_1 = 0.0f;
    float y_1 = 0.0f;
    float z_1 = 0.0f;

    float DispX = 0.0f;
    float DispY = 0.0f;
    float DispZ = 0.0f;

    float dt = delta_t;
    float Border;
    float Landscape;

    Quaternion draw_pos;

    sf::CircleShape shape;

public:

    ElementarElement(float x, float y, float z, float Vx, float Vy, float Vz, float landscape, float la);

    void Starting_conditions_for_W(float W[3]);

    void Move();

    void periodic_table();

    void Cout(int t);

    float Force(ElementarElement* element);

    void ThermoV();

    float KinEnergy();

    float* Imp(float p[3]);

    float Disp();

    float absV();

    

    void CoordFor3D(Quaternion camerapos, Quaternion camerarotation);

    void draw(sf::RenderWindow& window);

};
float moving(std::list<ElementarElement> ElementList, float Border, int t);
