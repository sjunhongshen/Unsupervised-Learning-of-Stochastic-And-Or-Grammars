#pragma once
#ifndef COMPONENTS_H
#define COMPONENTS_H

#include <unordered_map>
#include <string.h>
#include <iostream>
#include <opencv2/core/core.hpp>  
#include <opencv2/highgui/highgui.hpp>  
#include <opencv2/imgproc.hpp>
#include <random>
#include "face_params.h"

constexpr auto WINDOW_SIZE = 128;
constexpr auto PI = 3.1415926;

using namespace cv;
using namespace std;

class Components {
	/*
		x, y:	center of the components
		type:	component type (1, 2, 3...)
		scale:	a rescale coefficient
		orientation:	an angle it is rotated
	*/
	protected:
		int x = 0;
		int y = 0;
		int type = 1;
		int radius = 0;
		double scale = 1;
		double orientation = 0;
		bool RandFlag = true;
		bool negSampleFlag = true;

	public:
		Components(int type_in) { type = type_in; };

		vector<Point> getWindow();
		//	get/set methods
		void setX(int x);
		void setY(int y);
		void setType(int type);
		void setRadius(int radius);
		void setScale(double scale);
		void setOrientation(double orientation);
		void setRandFlag(bool Flag);
		void setNegSampleFlag(bool Flag);
		int getX();
		int getY();
		int getType();
		int getRadius();
		double getScale();
		double getOrientation();
		void copyXYSO(const vector<double> XYSO)
		{
			this->x = XYSO[0];
			this->y = XYSO[1];
			this->scale = XYSO[2];
			this->orientation = XYSO[3];
		};
		virtual Mat drawPic();
};

// Face
class Face : public Components {
	int w = WINDOW_SIZE;

	public:
		Face(int type) : Components(type) {
			this->x =  face_x_mean;
			this->y =  face_y_mean;
			this->radius = face_radius;
			this->negSampleFlag = false;
			// cout << "Face Initialized!\n"; 
		};
		Mat drawPic(Mat panel);
};

// Eye
class Eye : public Components {
	int w = WINDOW_SIZE;
	string eyeType;
	public:
		Eye(int type, string eyeType = "Left") : Components(type) {
			this->eyeType = eyeType;
			this->x =  eye_left_x_mean;
			this->y =  eye_left_y_mean;
			this->radius = eye_left_radius;
			this->negSampleFlag = false;
			if (this->eyeType == "Right"){
				this->x =  eye_right_x_mean;
				this->y =  eye_right_y_mean;
				this->radius = eye_right_radius;
			}
			// cout << eyeType <<" Eye Initialized!\n";
		};
		Mat drawPic(Mat panel);
};


// Ear
class Ear : public Components {
	int w = WINDOW_SIZE;
	string earType;
	public:
		Ear(int type, string earType = "Left") : Components(type) {
			this->earType = earType;
			this->x =  ear_left_x_mean;
			this->y =  ear_left_y_mean;
			this->radius = ear_left_radius;
			this->negSampleFlag = false;
			if (this->earType == "Right") {
				this->x =  ear_right_x_mean;
				this->y =  ear_right_y_mean;
				this->radius = ear_right_radius;
			}
			// cout << earType << " Ear Initialized!\n";
		};
		Mat drawPic(Mat panel);
};

// Nose
class Nose : public Components {
	int w = WINDOW_SIZE;
	public:
		Nose(int type) : Components(type) {
			this->x =  nose_x_mean;
			this->y =  nose_y_mean;
			this->radius = nose_radius;
			this->negSampleFlag = false;
			// cout << "Nose Initialized!\n";
		};
		Mat drawPic(Mat panel);
};

// Mouth
class Mouth : public Components {
	int w = WINDOW_SIZE;
	public:
		Mouth(int type) : Components(type) {
			this->x =  mouth_x_mean;
			this->y =  mouth_y_mean;
			this->radius = mouth_radius;
			this->negSampleFlag = false;
			// cout << "Mouth Initialized!\n";
		};
		Mat drawPic(Mat panel);
};


Point rescale(Point point, Point origin, double scale);
Point rotate(Point point, Point origin, double angle);
Point rescaleAndRotate(Point point, Point origin, double scale, double angle);
vector<Point> mirror(vector<Point> pointLists, int centerX);
Mat cropWindow(Mat panel, vector<Point> pointLists);

struct GeometricProps
{
	int x_, y_;
	double scale_, orientation_;
	double x_var_, y_var_;
	double scale_var_, orientation_var_;

	GeometricProps(int x, int y, double x_var, double y_var,
				   double scale_var, double orientation_var)
	:x_(x), y_(y), scale_(1.0), orientation_(0.0), x_var_(x_var),
	 y_var_(y_var), scale_var_(scale_var), orientation_var_(orientation_var)
	{}
};
#endif