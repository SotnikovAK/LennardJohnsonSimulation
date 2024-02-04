#pragma once
#include "libraries.h"
#include "Constants.h"

class Quaternion {
public:
	double w, x, y, z;

	Quaternion() {};

	Quaternion(double nx, double ny, double nz) {
		w = 0;
		x = nx;
		y = ny;
		z = nz;
	}

	Quaternion(double nw, double nx, double ny, double nz) {
		w = nw;
		x = nx;
		y = ny;
		z = nz;
	}

	Quaternion operator+(Quaternion b) {
		return Quaternion(w + b.w, x + b.x, y + b.y, z + b.z);
	}

	Quaternion operator-(Quaternion b) {
		return Quaternion(w - b.w, x - b.x, y - b.y, z - b.z);
	}

	Quaternion operator*(double b) {
		return Quaternion(b * w, b * x, b * y, b * z);
	}

	Quaternion operator*(Quaternion b) {

		return Quaternion(w * b.w - x * b.x - y * b.y - z * b.z,
			w * b.x + x * b.w + y * b.z - z * b.y,
			w * b.y - x * b.z + y * b.w + z * b.x,
			w * b.z + x * b.y - y * b.x + z * b.w);
	}

	Quaternion operator/(double b) {
		return Quaternion(w / b, x / b, y / b, z / b);
	}

	Quaternion inverse() {
		return Quaternion(w, -x, -y, -z) / get_magnitude();
	}


	double coord(int flag) {

		double A;

		switch (flag) {
		case 1:
			A = x;
			break;
		case 2:
			A = y;
			break;
		case 3:
			A = z;
			break;
		}

		return A;

	}
	void set_coord(float nx, float ny, float nz) {

		x = nx;
		y = ny;
		z = nz;

	}

	double get_magnitude() {
		return std::pow(std::pow(w, 2) + std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2), 0.5);
	}

	void set_magnitude(double newmag) {
		double factor = newmag / get_magnitude();

		w = w * factor;
		x = x * factor;
		y = y * factor;
		z = z * factor;

	}

	void normalize() {
		set_magnitude(1.);
	}

	Quaternion normalized() {
		Quaternion Q(w, x, y, z);
		Q.normalize();
		return Q;
	}

	sf::Vector2f getScreenPos() {
		if (z == 0)
			return sf::Vector2f(-1000, -1000);
		return sf::Vector2f(
			ratio * (screenwidth / 2.0) * x / z + (screenwidth / 2.0),
			ratio * (screenheight / 2.0) * y / z + (screenheight / 2.0));
	}
};