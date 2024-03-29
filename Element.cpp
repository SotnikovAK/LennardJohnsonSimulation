#include "Element.h"
#include "Quaternion.h"



ElementarElement::ElementarElement(float x_, float y_, float z_, float Vx_, float Vy_, float Vz_, float land_, float _la)
{
	std::tie(x, y, z) = std::make_tuple(x_, y_, z_);
	std::tie(Vx, Vy, Vz) = std::make_tuple(Vx_, Vy_, Vz_);

	std::tie(x_1, y_1, z_1) = std::make_tuple(x + Vx * dt, y + Vy * dt, z + Vz * dt);
	std::tie(DispX, DispY, DispZ) = std::make_tuple(0.0f, 0.0f, 0.0f);
	std::tie(Wx, Wy, Wz) = std::make_tuple(0.0f, 0.0f, 0.0f);

	dt = delta_t;
	Border = land_;
	Landscape = _la;

}

void ElementarElement::Starting_conditions_for_W(float W[3])
{
	Wx = W[0];
	Wy = W[1];
	Wz = W[2];
}


void ElementarElement::Move()
{

	float temp_cooard[3] = { x ,y ,z };
	float squared_time = pow(dt, 2);

	std::tie(x, y, z) = std::make_tuple(2 * x - x_1 + Wx * squared_time, 2 * y - y_1 + Wy * squared_time, 2 * z - z_1 + Wz * squared_time);
	std::tie(Vx, Vy, Vz) = std::make_tuple(0.5 * (x - x_1) / dt, 0.5 * (y - y_1) / dt, 0.5 * (z - z_1) / dt);

	std::tie(x_1, y_1, z_1) = std::make_tuple(temp_cooard[0], temp_cooard[1], temp_cooard[2]);

	std::tie(DispX, DispY, DispZ) = std::make_tuple(Vx, Vy, Vz);
}

void ElementarElement::periodic_table()
{
	float temp_coord[3] = { x ,y ,z };
	float temp_coord_[3] = { x_1 ,y_1 ,z_1 };

	for (int k = 0; k < 3; ++k) {
		if (temp_coord[k] < 0) {
			int temp = (-temp_coord[k] + Border - 1) / Border;
			temp_coord[k] += temp * Border;
			temp_coord_[k] += temp * Border;
		}
		else if (temp_coord[k] >= Border) {

			int temp = temp_coord[k] / Border;

			temp_coord[k] -= temp * Border;
			temp_coord_[k] -= temp * Border;
		}
	}


	std::tie(x, y, z) = std::make_tuple(temp_coord[0], temp_coord[1], temp_coord[2]);
	std::tie(x_1, y_1, z_1) = std::make_tuple(temp_coord_[0], temp_coord_[1], temp_coord_[2]);

}




void ElementarElement::Cout(int t)
{

	std::cout << std::endl << "time = " << t << '\n';
	std::cout << "x: " << x << " y: " << y << "z = " << z << '\n';
	std::cout << "|  Vx: " << Vx << " Vy : " << Vy << "Vz = " << Vz << '\n';
	std::cout << "| Wx : " << Wx << "Wy : " << Wy << "Wz = " << Wz << '\n';
}


void ElementarElement::ThermoV()
{

	std::seed_seq seed_{ time(NULL) };
	std::mt19937 mt(seed_);


	std::uniform_real_distribution<double> dis(0.0, 1);
	float w = dis(mt);

	float E_max = sqrt(bolzman * Temperature / MassParticale);
	std::normal_distribution<float> E(0, E_max);

	float k = sqrt(KT / MassParticale);

	if (w < Collisation) {
		std::tie(Vx, Vy, Vz) = std::make_tuple(E(mt) / k, E(mt) / k, E(mt) / k);
		std::tie(x_1, y_1, z_1) = std::make_tuple(x - Vx * dt, y - Vy * dt, z - Vz * dt);
	}

}

float ElementarElement::KinEnergy()
{
	float K = 0.0;
	K += (pow(Vx, 2) + pow(Vy, 2) + pow(Vz, 2)) / 2;
	return K;
}

float* ElementarElement::Imp(float p[3])
{
	std::tie(p[0], p[1], p[2]) = std::make_tuple(Vx, Vy, Vz);

	return p;
}

float ElementarElement::Disp()
{
	float D = 0.0;
	D += (pow(DispX, 2) + pow(DispY, 2) + pow(DispZ, 2)) / 2;
	return D;
}

float ElementarElement::absV()
{
	float vabs = 0.0;
	vabs += sqrt(pow(Vx, 2) + pow(Vy, 2) + pow(Vz, 2));
	return vabs;
}




void ElementarElement::CoordFor3D(Quaternion camerapos, Quaternion camerarotation)
{

	Quaternion pos(0, x, y, z);

	double distanceFromCamera;
	Quaternion draw_pos;

	draw_pos = camerarotation * ((pos - camerapos) * camerarotation.inverse());

	shape.setFillColor(sf::Color(std::rand() % 255, std::rand() % 255, std::rand() % 255));
	distanceFromCamera = draw_pos.get_magnitude();


	if (draw_pos.z < 0)
		shape.setRadius(0);
	else
		shape.setRadius(R / distanceFromCamera * 1000);

	shape.setPosition(draw_pos.getScreenPos() + sf::Vector2f(-R / distanceFromCamera * 1000 / Border, -R / distanceFromCamera * 1000 / Border));

}

void ElementarElement::draw(sf::RenderWindow& window)
{
	window.draw(shape);
}
