#pragma once

#include <string>

class Vector3 {
public:

	double x, y, z;

	Vector3(double x = 0.0, double y = 0.0, double z = 0.0)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}

	Vector3 operator+(const Vector3& v) {
		return Vector3(this->x + v.x, this->y + v.y, this->z + v.z);
	}

	Vector3 operator-(const Vector3& v) {
		return Vector3(this->x - v.x, this->y - v.y, this->z - v.z);
	}

	Vector3 operator+(const double c) {
		return Vector3(this->x + c, this->y + c, this->z + c);
	}

	Vector3 operator-(const double c) {
		return Vector3(this->x - c, this->y - c, this->z - c);
	}

	void operator+=(const Vector3& v) {
		this->x += v.x;
		this->y += v.y;
		this->z += v.z;
	}

	void operator-=(const Vector3& v) {
		this->x -= v.x;
		this->y -= v.y;
		this->z -= v.z;
	}

	void operator+=(const double c) {
		this->x += c;
		this->y += c;
		this->z += c;
	}

	void operator-=(const double c) {
		this->x -= c;
		this->y -= c;
		this->z -= c;
	}

	Vector3 operator*(const Vector3& v) {
		return Vector3(this->x * v.x, this->y * v.y, this->z * v.z);
	}

	Vector3 operator/(const Vector3& v) {
		return Vector3(this->x / v.x, this->y / v.y, this->z / v.z);
	}

	Vector3 operator*(const double c) {
		return Vector3(this->x * c, this->y * c, this->z * c);
	}

	Vector3 operator/(const double c) {
		return Vector3(this->x / c, this->y / c, this->z / c);
	}

	void operator*=(const Vector3& v) {
		this->x *= v.x;
		this->y *= v.y;
		this->z *= v.z;
	}

	void operator/=(const Vector3& v) {
		this->x /= v.x;
		this->y /= v.y;
		this->z /= v.z;
	}

	void operator*=(const double c) {
		this->x *= c;
		this->y *= c;
		this->z *= c;
	}

	void operator/=(const double c) {
		this->x /= c;
		this->y /= c;
		this->z /= c;
	}

	bool operator==(const Vector3& vec)
	{
		return this->x == vec.x && this->y == vec.y && this->z == vec.z;
	}

	bool operator!=(const Vector3& vec)
	{
		return !(this->x == vec.x && this->y == vec.y && this->z == vec.z);
	}

	std::string to_string() {
		return "(" + std::to_string(this->x) + ", " + std::to_string(this->y) + ", " + std::to_string(this->z) + ")";
	}
};


class Vector3i {
public:

	int x, y, z;

	Vector3i(int x = 0.0f, int y = 0.0f, int z = 0.0f)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}

	Vector3i(Vector3 vec)
	{
		this->x = (int)vec.x;
		this->y = (int)vec.y;
		this->z = (int)vec.z;
	}

	Vector3i operator+(const Vector3i& v) {
		return Vector3i(this->x + v.x, this->y + v.y, this->z + v.z);
	}

	Vector3i operator-(const Vector3i& v) {
		return Vector3i(this->x - v.x, this->y - v.y, this->z - v.z);
	}

	Vector3i operator+(const int c) {
		return Vector3i(this->x + c, this->y + c, this->z + c);
	}

	Vector3i operator-(const int c) {
		return Vector3i(this->x - c, this->y - c, this->z - c);
	}

	void operator+=(const Vector3i& v) {
		this->x += v.x;
		this->y += v.y;
		this->z += v.z;
	}

	void operator-=(const Vector3i& v) {
		this->x -= v.x;
		this->y -= v.y;
		this->z -= v.z;
	}

	void operator+=(const int c) {
		this->x += c;
		this->y += c;
		this->z += c;
	}

	void operator-=(const int c) {
		this->x -= c;
		this->y -= c;
		this->z -= c;
	}

	Vector3i operator*(const Vector3i& v) {
		return Vector3i(this->x * v.x, this->y * v.y, this->z * v.z);
	}

	Vector3i operator/(const Vector3i& v) {
		return Vector3i(this->x / v.x, this->y / v.y, this->z / v.z);
	}

	Vector3i operator*(const int c) {
		return Vector3i(this->x * c, this->y * c, this->z * c);
	}

	Vector3i operator/(const int c) {
		return Vector3i(this->x / c, this->y / c, this->z / c);
	}

	void operator*=(const Vector3i& v) {
		this->x *= v.x;
		this->y *= v.y;
		this->z *= v.z;
	}

	void operator/=(const Vector3i& v) {
		this->x /= v.x;
		this->y /= v.y;
		this->z /= v.z;
	}

	void operator*=(const int c) {
		this->x *= c;
		this->y *= c;
		this->z *= c;
	}

	void operator/=(const int c) {
		this->x /= c;
		this->y /= c;
		this->z /= c;
	}

	bool operator==(const Vector3i& vec)
	{
		return this->x == vec.x && this->y == vec.y && this->z == vec.z;
	}

	bool operator!=(const Vector3i& vec)
	{
		return !(this->x == vec.x && this->y == vec.y && this->z == vec.z);
	}

	std::string to_string() {
		return "(" + std::to_string(this->x) + ", " + std::to_string(this->y) + ", " + std::to_string(this->z) + ")";
	}
};
