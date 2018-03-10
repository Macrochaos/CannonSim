#pragma once
class Cannon
{
	struct Vec3 {
		union {
			float x, y, z = 0;
			float v[3];
		};
	};
	Vec3 Pos;
	Vec3 Speed;
	Vec3 Accel;
	struct Orientation{
		float theta;
		float phi;
	};
	
public:
	Orientation ThetaPhi;
	Cannon();
	~Cannon();
};

