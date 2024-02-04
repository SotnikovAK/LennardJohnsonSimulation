#pragma once
#pragma once

#include "libraries.h"
#include "Constants.h"
#include "Quaternion.h"

class ElementarElement
{
private:
    float M = MassParticale;

    float Vx = 0;
    float Vy = 0;
    float Vz = 0;

    float Wx = 0;
    float Wy = 0;
    float Wz = 0;

    float R = 0.1;

    float x = 0;
    float y = 0;
    float z = 0;

    float x_1 = 0;
    float y_1 = 0;
    float z_1= 0;

    float DispX, DispY, DispZ;

    float dt = delta_t;
    float Border;
    float Landscape;

    Quaternion draw_pos;

    sf::CircleShape shape;

public:

    ElementarElement(float x, float y, float z, float Vx, float Vy, float Vz, float landscape, float la);

    void Starting_conditions_for_W();

    void Move();

    void PeriodCoord();

    void Cout();

    float Force(ElementarElement* element);

    void ThermoV();

    float KinEnergy();

    float* Imp(float p[3]);

    float Disp();

    float absV();


    void CoordFor3D(Quaternion camerapos, Quaternion camerarotation);

    void draw(sf::RenderWindow& window);

};
