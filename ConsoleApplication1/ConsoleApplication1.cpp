// ConsoleApplication1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <vector>
#include <string>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include <math.h>
#include "Cannon.h"

using namespace std;
#define Rad_2_Deg 57.295779513
#define Deg_2_Rad 0.01745329252
#define PI 3.141592

#define QuadraticMinDist 40000
#define MinDist 200
#define FrameTime_Sec 0.01

#define Grav -9.8
#define WindSpeed -9.0
#define WindAcceleration -0.0995
#define Dist_T2T 5.0f

#define ShotSpeed 30
#define ShotSpeed_Needed 46.5
#define BallMass_kg 20
#define Drag_OnBall_sans_sqrVel 5.074
#define Blast_Radius 3.0

//
bool find_uint_InVector(vector<unsigned int> vec, int val) {
	int s = vec.size();
	for (int i = 0; i < s; ++i) {
		if (vec[i] = val)
			return true;
	}
	return false;
}
//------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------
//////////////////////////////////////////// v Struct Vec3 v ////////////////////////////////////////////
struct Vec3 {
	union {
		float x, y, z = 0;
		float v[3];
	};
	Vec3() {
		x = 0; y = 0; z = 0;
	}
	Vec3(float xx, float yy, float zz) {
		x = xx; y = yy; z = zz;
	}
	Vec3 operator = (Vec3 other) {
		x = other.x; y = other.y; z = other.z;
		return *this;
	}
	Vec3 operator +=(Vec3 other) {
		x += other.x;
		y += other.y;
		z += other.z;
		return *this;
	}
	Vec3 operator *(float f) {
		x *= f;
		y *= f;
		z *= f;
		return *this;
	}
	Vec3 operator *=(float f) {
		x *= f;
		y *= f;
		z *= f;
		return *this;
	}
	Vec3 operator /(float f) {
		x /= f;
		y /= f;
		z /= f;
		return *this;
	}
	Vec3 operator -(Vec3 other) {
		x -= other.x;
		y -= other.y;
		z -= other.z;
		return *this;
	}
	Vec3 operator +(Vec3 other) {
		x += other.x;
		y += other.y;
		z += other.z;
		return *this;
	}
	void Normalize(Vec3 &vect) {
		float mag = 0;
		mag = vect.x * vect.x + vect.y * vect.y + vect.z * vect.z;
		mag = sqrt(mag);
		vect.x *= mag; vect.y *= mag; vect.z *= mag;
	}
	float Mag() {
		return sqrt(x*x + y*y + z*z);
	}
	//float
};
Vec3 vecMidPoint(Vec3 &head, Vec3 &tail) {
	return (head - tail) * 0.5;
}
Vec3 Normalized(Vec3 &vec) {
	return vec *( 1 / vec.Mag());
}
float Dot(Vec3 & vecA,Vec3 & vecB) {
	Vec3 vecN_A = Normalized(vecA);
	Vec3 vecN_B = Normalized(vecB);
	return(vecN_A.x*vecN_B.x + vecN_A.y*vecN_B.y + vecN_A.z*vecN_B.z);
}

//////////////////////////////////////////// ^ Struct Vec3 ^ ////////////////////////////////////////////
//////////////////////////////////////////// v Struct Orientation v ////////////////////////////////////////////////
struct Orientation {
	float theta;
	float phi;
};
//////////////////////////////////////////// ^ Struct Orientation ^ ////////////////////////////////////////////////
//////////////////////////// v Struct Triangle v ////////////////////////////
struct Triangle {
	Vec3 V0;
	Vec3 V1;
	Vec3 V2;
	Vec3 cent;
	Triangle() {
		V0.x = V0.y = V0.z = 0;
		V1.x = Dist_T2T; V1.y = V1.z = 0;
		V2.x = Dist_T2T*0.5; V1.y = 0; V2.z = - (sqrt((Dist_T2T*Dist_T2T) - ((Dist_T2T*0.5)*((Dist_T2T*0.5)))));
		cent = vecMidPoint(V2, vecMidPoint(V1, V0));
	}
	Triangle operator=(Triangle other) {
		V0 = other.V0; V1 = other.V1; V2 = other.V2; cent = other.cent;
	}
	Triangle (Vec3 v_0) {
		V0 = v_0;
		V1.x = v_0.x + Dist_T2T; V1.y = v_0.y; V1.z = v_0.z;
		V2.x = v_0.x + Dist_T2T*0.5; V1.y = v_0.y; V2.z = v_0.z - (sqrt((Dist_T2T*Dist_T2T) - ((Dist_T2T*0.5)*((Dist_T2T*0.5)))));
		cent = (v_0 + vecMidPoint(V2, vecMidPoint(V1, V0)));
	}
	void flip_tri() {
		V1.z += 2 * (sqrt((Dist_T2T*Dist_T2T) - ((Dist_T2T*0.5)*((Dist_T2T*0.5)))));
	}
};

//////////////////////////// ^ Struct Triangle ^ ////////////////////////////
////////////////////////////////////// v Init funcs v //////////////////////////////////////
void init_0(Vec3 &v) {
	v.x = v.y = v.z = 0;
}
void init_Shoot(Vec3 &Speed, Cannon &boomy) {
	Speed.x = cos(boomy.ThetaPhi.phi) * ShotSpeed;
	Speed.z = sin(boomy.ThetaPhi.phi) * ShotSpeed;
	Speed.y = sin(boomy.ThetaPhi.theta) * ShotSpeed;
}
void init_firstT(Vec3 & first_target) {
	init_0(first_target);
	first_target.x = 200;
}
void init_vals(Vec3 &v, float x, float z) {
	init_0(v);
	v.x = x;
	v.z = z;
}
void init_targets20(Vec3 & first_target, Vec3* targetsVector) {
	for (int i = 0; i < 20; ++i) {
		init_0(targetsVector[i]);
	}
	init_firstT(targetsVector[0]);
	init_vals(targetsVector[9], 200, 4.33013 * 2);
	init_vals(targetsVector[18], 200, 4.33013 * 4);

	init_vals(targetsVector[1], 200 + (2.5 * 1), 4.33013 * 1);
	init_vals(targetsVector[10], 200 + (2.5 * 1), 4.33013 * 3);

	init_vals(targetsVector[2], 200 + (2.5 * 2), 4.33013 * 0);
	init_vals(targetsVector[11], 200 + (2.5 * 2), 4.33013 * 2);
	init_vals(targetsVector[19], 200 + (2.5 * 2), 4.33013 * 4);

	init_vals(targetsVector[3], 200 + (2.5 * 3), 4.33013 * 1);
	init_vals(targetsVector[12], 200 + (2.5 * 3), 4.33013 * 3);

	init_vals(targetsVector[4], 200 + (2.5 * 4), 4.33013 * 0);
	init_vals(targetsVector[13], 200 + (2.5 * 4), 4.33013 * 2);

	init_vals(targetsVector[5], 200 + (2.5 * 5), 4.33013 * 1);
	init_vals(targetsVector[14], 200 + (2.5 * 5), 4.33013 * 3);

	init_vals(targetsVector[6], 200 + (2.5 * 6), 4.33013 * 0);
	init_vals(targetsVector[15], 200 + (2.5 * 6), 4.33013 * 2);

	init_vals(targetsVector[7], 200 + (2.5 * 7), 4.33013 * 1);
	init_vals(targetsVector[16], 200 + (2.5 * 7), 4.33013 * 3);

	init_vals(targetsVector[8], 200 + (2.5 * 8), 4.33013 * 0);
	init_vals(targetsVector[17], 200 + (2.5 * 8), 4.33013 * 2);
}

void init_targets_tris(Vec3 & first_target, int &columns, int &targets, vector<Vec3> &Tvect, vector<Triangle> & TriVec) {
	int currRow = 0;
	int tris = targets / 3;
	int rows = tris / columns;
	float z_offset = (sqrt((Dist_T2T*Dist_T2T) - ((Dist_T2T*0.5)*((Dist_T2T*0.5)))));
	Vec3 aux;
	aux.z = z_offset * 2;
	//
	//Creates first triangle
	//
	Triangle trii(first_target);
	TriVec.push_back(trii);
	for (int k = 1; k < tris; ++k) {
		currRow = k / columns;
		if (k % 2 == 1) {
			
			Triangle trii(TriVec[k - 1].V1 + aux * ((k / columns) * 2));
			trii.flip_tri();
			TriVec.push_back(trii);
		}
		else {
			Triangle trii(TriVec[k - 1].V1 + aux * ((k / columns) * 2));
			TriVec.push_back(trii);
		}
	}
	//put all Vec3 from triangles in a vector
	for (int i = 0; i < TriVec.size(); ++i) {
		Tvect.push_back(TriVec[i].V0);
		Tvect.push_back(TriVec[i].V1);
		Tvect.push_back(TriVec[i].V2);
	}
	int remainingTargs = targets - TriVec.size() * 3;
	if (remainingTargs == 1) {
		Vec3 remVec1 = first_target + (aux*currRow);
		Tvect.push_back(remVec1);
	}
	if (remainingTargs == 2) {
		Vec3 remVec1 = first_target + (aux*(currRow + 1));
		Tvect.push_back(remVec1);
		Vec3 remVec2 = first_target + (aux*(currRow + 1));
		remVec2.x += Dist_T2T;
		Tvect.push_back(remVec2);
	}
}


////////////////////////////////////// ^ Init funcs ^ //////////////////////////////////////
////////////////////////////////////// v Aim funcs v //////////////////////////////////////
float calcShot_time(const float &theta_ang) {
	return ((2 * ShotSpeed*sin(theta_ang)) / Grav);
}
float clacQuadrDist(Vec3 v) {
	return v.x * v.x + v.z * v.z;
}
float calcDistance(Vec3 v) {
	return sqrt(clacQuadrDist(v));
}
bool checkTripleKill(Vec3 v1, Vec3 v2, Vec3 v3) {
	float d1 = calcDistance(v1);
	float d2 = calcDistance(v2);
	float d3 = calcDistance(v3);
	if ((d2 - d1) <= 5.01 && (d3 - d1) <= 5.01 && (d3 - d2) <= 5.01)
		return true;
	return false;
}
void AimShot20(Vec3* vecvec, Cannon cann, float Time) {
	vector<unsigned int> usedTargets;
	usedTargets.push_back(0);
	Vec3 MemP1;				   float P1QD = 40000;
	Vec3 MemP2; init_0(MemP2); float P2QD = 0;
	Vec3 MemP3; init_0(MemP3); float P3QD = 0;

	Vec3 target;

	//Gets 3 closest points to
	MemP1 = vecvec[0];
	MemP2 = vecvec[1];
	for (int i = 2; i < 20; ++i) {
		if (clacQuadrDist(vecvec[i]) < clacQuadrDist(vecvec[1]))
			MemP2 = vecvec[i]; usedTargets.push_back(i);
	}
	for (int i = 3; i < 20; ++i) {
		if (clacQuadrDist(vecvec[i]) < clacQuadrDist(vecvec[2]) || !find_uint_InVector(usedTargets,i))
			MemP3 = vecvec[i]; usedTargets.push_back(i);
	}
	if (checkTripleKill(MemP1, MemP2, MemP3))
		target = vecMidPoint(MemP3, vecMidPoint(MemP2, MemP1));

}
float calcTheta(float & Mag) {
	return (asin(((Grav*Mag) / (ShotSpeed_Needed*ShotSpeed_Needed))) / 2);
}
void calcDir(Vec3 vec, Cannon cann, float Time) {
	float mag = vec.Mag();
	float currPhi = cann.ThetaPhi.phi;
	float currTheta = cann.ThetaPhi.theta;
	Vec3 Dir = Normalized(vec);
	float NTheta = calcTheta(mag);
	//
	// Calculate z correctinf due to wind
	//
	Vec3 WindCorrectedVec = vec;
	float flyTime = calcShot_time(NTheta);
	WindCorrectedVec.z += flyTime*-WindAcceleration;
	Vec3 CorrectionVec = WindCorrectedVec - vec;
	CorrectionVec = CorrectionVec * Dot(WindCorrectedVec, vec);
	float NMag = ((vec+ CorrectionVec).Mag());
	Dir = Normalized(vec + CorrectionVec);
	NTheta = calcTheta(NMag);
	flyTime = calcShot_time(NTheta);
	//
	//
	//
	float NPhi = atan(Dir.z / Dir.x);
	Time += ((NPhi - currPhi)*Rad_2_Deg) + ((NTheta - currTheta)*Rad_2_Deg) + flyTime;
	cann.ThetaPhi.theta = NTheta;
	cann.ThetaPhi.phi = NPhi;
}
void AimShot(vector<Vec3> vecvec, vector<Triangle> trivec, Cannon cann, float Time) {
	vector<Vec3> usedvecs;
	for (int i = 0; i < trivec.size(); ++i) {
		calcDir(trivec[i].cent, cann, Time);
		usedvecs.push_back(trivec[i].V0);
		usedvecs.push_back(trivec[i].V1);
		usedvecs.push_back(trivec[i].V2);
	}
	for (int i = usedvecs.size() - 1; i < vecvec.size(); ++i) {
		calcDir(vecvec[i], cann, Time);
	}
}
Vec3 calcShot_pos(const Orientation &aim, float Time) {
	Vec3 pointHit;
	Vec3 Dir;
	//basic aiming v
	Dir.x = cos(aim.phi);
	Dir.z = sin(aim.phi);
	float airTime = calcShot_time(aim.theta);
	Time += airTime;
	//basic aiming ^

	pointHit.x = airTime*Dir.x * ShotSpeed;
	pointHit.z = airTime*Dir.z * ShotSpeed;
	//add z displacement due to wind


	return pointHit;
}
////////////////////////////////////// ^ Aim funcs ^ //////////////////////////////////////
////////////////////////////////////// v Wind funcs v //////////////////////////////////////
void apply_Wind_Force(Vec3 Speed) {
	Speed.z += WindAcceleration * FrameTime_Sec;
}
////////////////////////////////////// ^ Wind funcs ^ //////////////////////////////////////
////////////////////////////////////// v Shot funcs v //////////////////////////////////////
bool shotHit(Vec3 &Pos) {
	if (Pos.y <= 0)
		return true;
	else
		return false;
}
void ShotFlies(Vec3 SpeedVec, Vec3 PosVec, float Time) {
	while (PosVec.y > 0) {
		PosVec += SpeedVec * FrameTime_Sec;
		apply_Wind_Force(SpeedVec);
		Time += FrameTime_Sec;
	}
}
////////////////////////////////////// ^ Shot funcs ^ //////////////////////////////////////
////////////////////////////////////// v Targets funcs v //////////////////////////////////////
bool is_quadr_dist(Vec3 & point) {
	if ((point.x * point.x + point.y* point.y) >= QuadraticMinDist)
		return true;
	else
		return false;
}

////////////////////////////////////// ^ Targets funcs ^ //////////////////////////////////////
//------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------
/////////////////////////////////////////// v Vital Variables v /////////////////////////////////////////////
float Time = 0.0000f;
int cols = 3;
int targets = 20;
vector<Vec3> targetsVector;
vector<Triangle> trivec;
Cannon Boom;
Vec3 firstT;
/////////////////////////////////////////// ^ Vital Variables ^ /////////////////////////////////////////////
//------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------
using namespace std;
/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

//func toma lista de ints los suma y regresa res for, while, recursivo

//func hace division recursiva sin usar '/'

//func que cheque si un string es un palindromo

//func tome 2 lista de string e intercale los chars

//func recibe lista de ints y regresa en entero más grande que se pueda generar con ellos

////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/
vector<int> v_int = { 1,3,5,7,9,11,13,15 };
vector<int> v_int2 = { 21,34,5,7,59,11,613,15 };
vector<int> v_int3 = { 1,3,25,7,9,121,13,155 };

int func1f(vector<int> v_int) {
	int res= 0;
	int vSize = v_int.size();
	for (int i = 0; i < vSize; ++i) {
		res += v_int[i];
	}
	return res;
}
int func1w(vector<int> v_int) {
	int sum = 0;
	//int aux=0;
	int s = v_int.size();
	while ((--s)>=0) {
		sum += v_int[s];
		//++aux;
	}
	return sum;
}
int func1r(vector<int> v_int) {
	static int sum = 0;
	static int aux = 0;
	static int vSize = v_int.size();
	if (aux < vSize) {
		sum += v_int[aux];
		++aux;
		func1r(v_int);
	}
	return sum;
}
int iSum(const vector<int>& vals) {
	if (vals.size()) {
		int res = vals[0];
		//int vector<int>newVals(vals[1]);
		if (vals.size() > 1) {
			//vector<int> nVec(vals[1], vals.size()-1);
			vector<int> nVec;
			nVec.insert(nVec.end(), vals.begin() + 1, vals.end());
			res += iSum(nVec);
		}
		return res;
	}
	//int res = vals[0];

	//res += iSum(vector<int> v_i = {2,3,4,5});
	return 0;
}
//////////////////////////////////////////////////////
float div_wo_div_r(int a, int b) {
	static float res = 0;
	static int a_mod = a;
	static int aux = 0;
	if (a_mod > b) {
		res += a - b;
		a_mod -= b;
		div_wo_div_r(a, b);
	}
	else if (a_mod < b) {
		++aux;
		res += b*pow(10, -aux);
		a_mod -= b*pow(10, -aux);
		div_wo_div_r(a, b);
	}
	return res;
}
int div_wo_div_i(int a, int b) {
	int res = 0;
	while (a > b) {
		++res;
		a -= b;
	}
	return res;
}

int div_wo_div_fast(int a, int b) {

	string num = to_string(b);
	string res = "";
	const int ascii_to_int = -48;
	int div_res = 0;
	int l = num.length();
	int digit = 0;
	int a_copy = a;
	int b_copy = b;
	
	for (int i = 0; i < l; ++i) {
		digit = static_cast<int>(num[i] - ascii_to_int);

	}

	/*int res = 0;
	while (a > b) {
		++res;
		a -= b;
	}*/
	return div_res;
}

int divdiv(int divid, int divis) {
	int res = 0;
	if (divis == 0) {
		assert(divis,"division by 0"); 
	}
	if (divid >= divis) {
		divid -= divis;
		res += divdiv(divid, divis) + 1;
	}
	return res;
}
//////////////////////////////////////////////////////
void reverse_STR(string &str) {
	int l = str.length();
	for (int i = 0; i < l / 2; ++i) {
		swap(str[i], str[l - i - 1]);
	}
}

bool is_palindr(string &s) {
	s.erase(remove(s.begin(), s.end(), ' '),s.end());
	int l = s.size();
	for (int i = 0; i < l; ++i) {
		s[i] = tolower(s[i]);
	}
	string s2 = s;
	reverse_STR(s2);
	if (s2==s)
		return true;
	else
		return false;
}

vector<int> int_number_to_digits(const int &num) {
	vector<int> ret_val;
	return ret_val;
}
//////////////////////////////////////////////////////
string str_mixer(string a, string b) {
	string c;
	int s1 = sizeof(a);
	int s2 = sizeof(b);
	if (a < b) {
		for (int i = 0; i < s1; ++i) {
			c[(i * 2) + 1] = a[i];
			c[i * 2] = b[i];
		}
		for (int i = s1+1; i < s2; ++i) {
			c[i] = b[i];
		}
	}
	else {
		for (int i = 0; i < s2; ++i) {
			c[(i * 2) + 1] = b[i];
			c[i * 2] = a[i];
		}
		for (int i = s2 + 1; i < s1; ++i) {
			c[i] = a[i];
		}
	}
	return c;
}
//////////////////////////////////////////////////////
unsigned int highest_int(vector<int> v_int) {
	//modules & divisions
	return 0;
}
//////////////////////////////////////////////////////
int Ackerman(unsigned long long m, unsigned long long n)
{
	if (m == 0)
		return n + 1;
	else
	{
		if (n == 0)
			return Ackerman(m - 1, 1);
		else
			return Ackerman(m - 1, Ackerman(m, n - 1));
	}
}

int main()
{
	/*
	unsigned long long m =			3		;
	unsigned long long n =			9		;
	cout << func1r(v_int2) << endl;
	cout << func1r(v_int) << endl;
	cout << func1w(v_int2) << endl;
	cout << func1f(v_int2) << endl;
	cout << iSum(v_int2) << endl;

	cout << divdiv(21, 7) << endl;
	string s = "  lo Ongg nooL";
	if (is_palindr(s))
		cout << "String is palindrome" << endl;
	else
		cout << "Error" << endl;

	cout << to_string(1234) << endl;

	cout <<"Ackerman number w/ M: "<<m<<" N: "<< n <<"	::	"<< Ackerman(m, n) << endl;

	float Time = 0.00f;
	Vec3 Wind;
	Vec3 BallSpeedVector;
	Vec3 BallPosVector;
	vector<Vec3> targetsVector;
	Cannon Boom;
	Vec3 firstT;
	*/
	////////////////// Init section //////////////////
	init_firstT(firstT);
	init_targets_tris(firstT, cols, targets, targetsVector,trivec);
	AimShot(targetsVector, trivec, Boom, Time);
	cout << "Time elapsed: " << Time << endl;

	////////////////// Init section //////////////////
    return 0;
}

