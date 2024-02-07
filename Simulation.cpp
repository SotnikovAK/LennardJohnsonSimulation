#include "Simulation.h"

#include "Simulation.h"
#include "Quaternion.h"

const sf::Vector2i size(screenwidth, screenheight);

Simulation::Simulation(float BorderSize, int n, int flag_animation_)
{
    CreationElements(BorderSize, n);
    flag_animation = flag_animation_;
}

static float ForceDependence(float ro) {
    float f1 = 0.0;
    float f2 = 0.0;

    f1 = 48 * pow(alpha, 12) / pow(ro, 7);
    f2 = -24 * pow(alpha, 6) / pow(ro, 4);

    return f1 + f2;

}
static float force(float delta_coord[3]) {
    float absF = 0.0;
    if (pow(delta_coord[0], 2) + pow(delta_coord[1], 2) + pow(delta_coord[2], 2) < FinishForceRo and pow(delta_coord[0], 2) + pow(delta_coord[1], 2) + pow(delta_coord[2], 2) > 0) {
        absF = ForceDependence(pow(delta_coord[0], 2) + pow(delta_coord[1], 2) + pow(delta_coord[2], 2));
    }
    return absF;
}



void Simulation::WindowInteraction()
{
    if (flag_animation) {

        sf::RenderWindow window(sf::VideoMode(screenwidth, screenheight), "Simulation");
        sf::Texture T;
        sf::Event event;

        SpecialWritingFunction s;

        T.create(window.getSize().x, window.getSize().y);

        window.setMouseCursorVisible(false);

        Quaternion camerapos(0, 0, 0, 0);
        Quaternion cameravel(0, 0, 0, 0);
        Quaternion camerarotation(Quaternion(1, 0, 0, 0).normalized());
        Quaternion temprotation(0, 0, 0, 0);

        int oldMouseX = screenwidth / 2;
        int oldMouseY = screenheight / 2;
        int dMouseX;
        int dMouseY;

        long long int  timer = 0;
        bool flag = true;

        while (flag)
        {
            if (!window.isOpen()) {
                flag = false;
            }

            while (window.pollEvent(event))
            {
                if (event.type == sf::Event::Closed or sf::Keyboard::isKeyPressed(sf::Keyboard::Escape))
                {
                    window.close();
                }

                if (event.type == sf::Event::MouseMoved) {
                    dMouseX = event.mouseMove.x - oldMouseX;
                    dMouseY = event.mouseMove.y - oldMouseY;

                    if ((dMouseX != 0) or (dMouseY != 0)) {

                        temprotation = Quaternion(1, 0.001 * dMouseY, -0.001 * dMouseX, 0).normalized();
                        camerarotation = temprotation * camerarotation;

                        sf::Mouse::setPosition(sf::Vector2i(screenwidth / 2, screenheight / 2), window);

                        oldMouseX = screenwidth / 2;
                        oldMouseY = screenheight / 2;
                    }


                }
            }

            /*
            for (auto element = ElementList.begin(); element != ElementList.end(); element++)
            {
                float W[3] = { 0.0,0.0,0.0 };
                element->Starting_conditions_for_W(W);
            }
            */
            U = 0.0;
            window.clear();
            float A[9][2] = {
                {-1.0, -1.0},
                {-1.0, 0.0},
                {-1.0, 1.0},
                {0.0, -1.0},
                {0.0,0.0},
                {0.0, 1.0},
                {1.0, -1.0},
                {1.0, 0.0},
                {1.0, 1.0},
            };

            float B[3] = { -1.0,0.0,1.0 };
            float absF = 0.0f;

            float X[max_n];
            float Y[max_n];
            float Z[max_n];

            int t = 0;
            for (auto element = ElementList.begin(); element != ElementList.end(); ++element) {
                for (int i = 0; i < 3; ++i) {
                    for (int k = 0; k < 9; ++k) {
                        std::tie(X[t], Y[t], Z[t]) = std::make_tuple(A[k][0] * BorderSize + element->x, A[k][1] * BorderSize + element->y, B[i] * BorderSize + element->z);
                        ++t;
                    }
                }
            }

            for (auto element = ElementList.begin(); element != ElementList.end(); ++element) {
                float W[3] = { 0.0f,0.0f,0.0f };
                for (int k = 0; k < t; ++k) {
                    float delta_coord[3] = { element->x - X[k], element->y - Y[k], element->z - Z[k] };
                    absF = force(delta_coord);
                    W[0] += absF * delta_coord[0];
                    W[1] += absF * delta_coord[1];
                    W[2] += absF * delta_coord[2];

                    if (pow(delta_coord[0], 2) + pow(delta_coord[1], 2) + pow(delta_coord[2], 2) > 0) {
                        U += 4 * (pow(alpha, 12) / pow(pow(delta_coord[0], 2) + pow(delta_coord[1], 2) + pow(delta_coord[2], 2), 6) - pow(alpha, 6) / pow(pow(delta_coord[0], 2) + pow(delta_coord[1], 2) + pow(delta_coord[2], 2), 3));
                    }

                }
                element->Starting_conditions_for_W(W);

            }

            for (auto element = ElementList.begin(); element != ElementList.end(); element++)
            {
                element->Cout(timer);
                element->Move();
                element->periodic_table();




                if (timer * delta_t < relaxTime * endtime) {

                    element->ThermoV();
                }

            }

            camerarotation.normalize();

            cameravel = cameravel + camerarotation.inverse() * (Quaternion(0, 0, 0, 0.001 * (double)sf::Keyboard::isKeyPressed(sf::Keyboard::W)) * camerarotation);
            cameravel = cameravel + camerarotation.inverse() * (Quaternion(0, 0, 0, -0.001 * (double)sf::Keyboard::isKeyPressed(sf::Keyboard::S)) * camerarotation);

            camerapos = camerapos + cameravel;

            if (sf::Keyboard::isKeyPressed(sf::Keyboard::Q)) {
                camerarotation = Quaternion(1, 0, 0, 0.01).normalized() * camerarotation;
            }
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::E)) {
                camerarotation = Quaternion(1, 0, 0, -0.01).normalized() * camerarotation;
            }




            for (auto element = ElementList.begin(); element != ElementList.end(); element++)
            {
                element->CoordFor3D(camerapos, camerarotation);
            }

            for (auto element = ElementList.begin(); element != ElementList.end(); element++)
            {
                element->draw(window);
            }
            if (timer % 16 == 0)
            {
                T.update(window);
            }


            window.display();

            ++timer;
        }
    }
    else {
        long long int  timer = 0;
        bool flag = true;
        float E[max_n][3];
        float V_[max_n];
        float D[max_n][2];
        float U = 0.0;
        float SumDisp = 0.0;
        SpecialWritingFunction s;
        while (flag)
        {
            if (timer > endtime) {
                flag = false;
                std::ofstream outputFile("E_.txt");
                if (outputFile.is_open()) {
                    for (int i = relaxTime + 1; i < timer; i++) {
                        for (int j = 0; j < 3; j++) {
                            outputFile << E[i][j] << " ";
                        }
                        outputFile << std::endl;
                    }
                    outputFile.close();
                    std::cout << "Содержимое массива успешно записано в файл." << std::endl;
                }

                std::ofstream outputFile2("V_.txt");
                if (outputFile2.is_open()) {
                    for (int i = relaxTime + 1; i < timer; i++) {
                         outputFile2 << V_[i] << " ";
                         outputFile2 << std::endl;
                    }
                    outputFile2.close();
                    std::cout << "Содержимое массива успешно записано в файл." << std::endl;
                }
                std::ofstream outputFile3("D_.txt");
                if (outputFile3.is_open()) {
                    for (int i = relaxTime + 1; i < timer; i++) {
                        for (int j = 0; j < 2; j++) {
                            outputFile3 << D[i][j] << " ";
                        }
                        outputFile3 << std::endl;
                    }
                    outputFile3.close();
                    std::cout << "Содержимое массива успешно записано в файл." << std::endl;
                }

            }

            U = 0.0;
            window.clear();
            float A[9][2] = {
                {-1.0, -1.0},
                {-1.0, 0.0},
                {-1.0, 1.0},
                {0.0, -1.0},
                {0.0,0.0},
                {0.0, 1.0},
                {1.0, -1.0},
                {1.0, 0.0},
                {1.0, 1.0},
            };

            float B[3] = { -1.0,0.0,1.0 };
            float absF = 0.0f;

            float X[max_n];
            float Y[max_n];
            float Z[max_n];

            int t = 0;
            for (auto element = ElementList.begin(); element != ElementList.end(); ++element) {
                for (int i = 0; i < 3; ++i) {
                    for (int k = 0; k < 9; ++k) {
                        std::tie(X[t], Y[t], Z[t]) = std::make_tuple(A[k][0] * BorderSize + element->x, A[k][1] * BorderSize + element->y, B[i] * BorderSize + element->z);
                        ++t;
                    }
                }
            }

            for (auto element = ElementList.begin(); element != ElementList.end(); ++element) {
                float W[3] = { 0.0f,0.0f,0.0f };
                for (int k = 0; k < t; ++k) {
                    float delta_coord[3] = { element->x - X[k], element->y - Y[k], element->z - Z[k] };
                    absF = force(delta_coord);
                    W[0] += absF * delta_coord[0];
                    W[1] += absF * delta_coord[1];
                    W[2] += absF * delta_coord[2];

                    if (pow(delta_coord[0], 2) + pow(delta_coord[1], 2) + pow(delta_coord[2], 2) > 0) {
                        U += 4 * (pow(alpha, 12) / pow(pow(delta_coord[0], 2) + pow(delta_coord[1], 2) + pow(delta_coord[2], 2), 6) - pow(alpha, 6) / pow(pow(delta_coord[0], 2) + pow(delta_coord[1], 2) + pow(delta_coord[2], 2), 3));
                    }


                }
                element->Starting_conditions_for_W(W);

            }

            for (auto element = ElementList.begin(); element != ElementList.end(); element++)
            {
                /*element->Cout(timer);*/
                element->Move();
                element->periodic_table();




                if (timer * delta_t < relaxTime * endtime) {

                    element->ThermoV();
                }
                if (timer > relaxTime) {
                    float SumKinEnergy = 0.0;
                    for (auto element = ElementList.begin(); element != ElementList.end(); element++) {
                        SumKinEnergy += element->KinEnergy();
                    }
                    E[timer][0] = timer;
                    E[timer][1] = SumKinEnergy;
                    E[timer][2] = U;
                    
                    for (auto element = ElementList.begin(); element != ElementList.end(); element++) {
                        V_[timer] = element->absV();
                    }
                }
                for (auto element = ElementList.begin(); element != ElementList.end(); element++) {
                    SumDisp += element->Disp();
                }
                D[timer][0] = timer;
                D[timer][1] = SumDisp;


            }

            ++timer;

        }
    }
}

static float random_zero_one() {
    return (double)rand() / RAND_MAX;
}
static float sign(float x) {
    return x / abs(x);
}

float MinDistance(float minro, float BorderSize, float delta_coord[3]) {
    for (int j = 0; j < 3; ++j) {
        if (BorderSize / 2 < abs(delta_coord[j])) {
            delta_coord[j] = sign(delta_coord[j]) * (BorderSize - delta_coord[j]);
        }

    }
    if (minro > pow(delta_coord[0], 2) + pow(delta_coord[1], 2) + pow(delta_coord[2], 2)) {
        minro = pow(delta_coord[0], 2) + pow(delta_coord[1], 2) + pow(delta_coord[2], 2);
    }

    return minro;

}

void Simulation::CreationElements(float BorderSize_, int N)
{
    BorderSize = BorderSize_;

    std::seed_seq seed_{ time(NULL) };
    std::mt19937 mt(seed_);
    landscape = screenwidth / BorderSize;

    std::cout << "W|S - camera_movement, Q|R - rotation of camera, ESC - exit" << std::endl;


    float minro = 10 * alpha;
    float r = 1.4 * alpha;


    float X[max_n];
    float Y[max_n];
    float Z[max_n];

    float Vx[max_n];
    float Vy[max_n];
    float Vz[max_n];


    for (int k = 0; k < N; ++k) {

        std::tie(X[k], Y[k], Z[k]) = std::make_tuple(BorderSize * random_zero_one(), BorderSize * random_zero_one(), BorderSize * random_zero_one());
        std::tie(Vx[k], Vy[k], Vz[k]) = std::make_tuple(2 * random_zero_one() - 1, 2 * random_zero_one() - 1, 2 * random_zero_one() - 1);

        for (int u = 0; u < k; ++u) {
            float delta_coord[3] = { X[k] - X[u] , Y[k] - Y[u] , Z[k] - Z[u] };

            minro = MinDistance(minro, BorderSize, delta_coord);

        }

        while (minro < r) {
            std::tie(X[k], Y[k], Z[k]) = std::make_tuple(BorderSize * random_zero_one(), BorderSize * random_zero_one(), BorderSize * random_zero_one());
            minro = 10 * alpha;

            for (int u = 0; u < k; ++u) {
                float delta_coord[3] = { X[k] - X[u] , Y[k] - Y[u] , Z[k] - Z[u] };

                minro = MinDistance(minro, BorderSize, delta_coord);
            }

        }




    }
    for (int i = 0; i < N; ++i) {
        ElementList.push_back(ElementarElement(X[i], Y[i], Z[i], Vx[i], Vy[i], Vz[i], BorderSize, landscape));

    }

}
