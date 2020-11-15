#include <iostream>  
#include <unordered_map>
#include <string.h>
#include <vector>
#include <opencv2/core/core.hpp>  
#include <opencv2/highgui/highgui.hpp>  
#include <opencv2/imgproc.hpp>

#include "components.h"
// Position variance
double face_x_var = 5;
double face_y_var = 5;
double nose_x_var = 3;
double nose_y_var = 3;
double mouth_x_var = 3;
double mouth_y_var = 3;
double ear_left_x_var = 3;
double ear_left_y_var = 3;
double ear_right_x_var = 3;
double ear_right_y_var = 3;
double eye_left_x_var = 3;
double eye_left_y_var = 3;
double eye_right_x_var = 3;
double eye_right_y_var = 3;

// Scale variance
double face_scale_var = 0.1;
double nose_scale_var = 0.1;
double mouth_scale_var = 0.1;
double ear_left_scale_var = 0.1;
double ear_right_scale_var = 0.1;
double eye_left_scale_var = 0.1;
double eye_right_scale_var = 0.1;

// Orientation variance
double face_pose_var = 10;
double nose_pose_var = 10;
double mouth_pose_var = 10;
double ear_left_pose_var = 10;
double ear_right_pose_var = 10;
double eye_left_pose_var = 10;
double eye_right_pose_var = 10;

// radius
int face_radius = 55;
int eye_left_radius = 12;
int eye_right_radius = 12;
int ear_left_radius = 12;
int ear_right_radius = 12;
int nose_radius = 17;
int mouth_radius = 17;


int thickness = 2;

// Negtive Samples range
int lower_limit = 20;
int higher_limit = 100;

std::default_random_engine generator;

std::normal_distribution<double> face_x_dist(face_x_mean, face_x_var);
std::normal_distribution<double> face_y_dist(face_y_mean, face_y_var);
std::normal_distribution<double> face_scale_dist(1.0, face_scale_var);
std::normal_distribution<double> face_orientation_dist(0.0, face_pose_var);

std::normal_distribution<double> nose_x_dist(nose_x_mean, nose_x_var);
std::normal_distribution<double> nose_y_dist(nose_y_mean, nose_y_var);
std::normal_distribution<double> nose_scale_dist(1.0, nose_scale_var);
std::normal_distribution<double> nose_orientation_dist(0.0, nose_pose_var);

std::normal_distribution<double> mouth_x_dist(mouth_x_mean, mouth_x_var);
std::normal_distribution<double> mouth_y_dist(mouth_y_mean, mouth_y_var);
std::normal_distribution<double> mouth_scale_dist(1.0, mouth_scale_var);
std::normal_distribution<double> mouth_orientation_dist(0.0, mouth_pose_var);

std::normal_distribution<double> ear_left_x_dist(ear_left_x_mean, ear_left_x_var);
std::normal_distribution<double> ear_left_y_dist(ear_left_y_mean, ear_left_y_var);
std::normal_distribution<double> ear_left_scale_dist(1.0, ear_left_scale_var);
std::normal_distribution<double> ear_left_orientation_dist(0.0, ear_left_pose_var);

std::normal_distribution<double> ear_right_x_dist(ear_right_x_mean, ear_right_x_var);
std::normal_distribution<double> ear_right_y_dist(ear_right_y_mean, ear_right_y_var);
std::normal_distribution<double> ear_right_scale_dist(1.0, ear_right_scale_var);
std::normal_distribution<double> ear_right_orientation_dist(0.0, ear_right_pose_var);

std::normal_distribution<double> eye_left_x_dist(eye_left_x_mean, eye_left_x_var);
std::normal_distribution<double> eye_left_y_dist(eye_left_y_mean, eye_left_y_var);
std::normal_distribution<double> eye_left_scale_dist(1.0, eye_left_scale_var);
std::normal_distribution<double> eye_left_orientation_dist(0.0, eye_left_pose_var);


std::normal_distribution<double> eye_right_x_dist(eye_right_x_mean, eye_right_x_var);
std::normal_distribution<double> eye_right_y_dist(eye_right_y_mean, eye_right_y_var);
std::normal_distribution<double> eye_right_scale_dist(1.0, eye_right_scale_var);
std::normal_distribution<double> eye_right_orientation_dist(0.0, eye_right_pose_var);

std::uniform_int_distribution<int> negSamples_x_dist(lower_limit, higher_limit);
std::uniform_int_distribution<int> negSamples_y_dist(lower_limit, higher_limit);

/*
	Inputs:
		point: point that needed to be rescaled
		orginal: the original point
		scale:	scale coefficient;
	Return:
		a point rescaled.
*/

Point rescale(Point point, Point origin, double scale) {
	Point diff = Point(point.x - origin.x, point.y - origin.y);
	double diff_x = diff.x * scale;
	double diff_y = diff.y * scale;
	return Point(round(origin.x + diff_x), round(origin.y + diff_y));
}

/*
	Inputs:
		point: point that needed to be rotated
		orginal: the original point
		angle:	the angle rotated (clock-wise)
	Return:
		a point rotated
*/


Point rotate(Point point, Point origin, double angle) {
	int x = point.x - origin.x;
	int y = point.y - origin.y;
	double distance = sqrt(x * x + y * y);
	double theta = atan2(y, x) * 180 / PI;
	double theta_new = theta + angle;

	double new_x = distance * cos(theta_new * PI / 180);
	double new_y = distance * sin(theta_new * PI / 180);
	return Point(round(origin.x + new_x), round(origin.y + new_y));
}

Point rescaleAndRotate(Point point, Point origin, double scale, double angle) {
	Point A = rescale(point, origin, scale);
	Point B = rotate(A, origin, angle);
	return B;
}


vector<Point> mirror(vector<Point> pointLists, int centerX) {
	vector<Point> newPointLists;
	for (int i = 0; i < pointLists.size(); i++) {
		Point A = pointLists[i];
		int x = centerX * 2 - A.x;
		newPointLists.push_back(Point(x, A.y));
	}
	return newPointLists;
}

/*
	Crop a window picture on a panel
	Input: vector<Point>, panel
		   4 points referring to the left-upper, right-upper, left-lower, right-lower point of a rectangle window
	output: Mat
			a window image cropped
*/

Mat cropWindow(Mat panel, vector<Point> pointLists) {
	int width = pointLists[1].x - pointLists[0].x;
	int height = pointLists[2].y - pointLists[0].y;
	Rect rect(pointLists[0].x, pointLists[0].y, width, height);
	return panel(rect);
}


/*
	get the window for a component
	Input: null
	output: a point list refering to the four points left-upper, right-upper, left-lower, right-lower
*/

vector<Point> Components::getWindow() {
	vector<Point> windowPoints;
	windowPoints.push_back(Point(round(x - radius * scale > 0 ? x - radius * scale : 0), round(y - radius * scale > 0 ? y - radius * scale : 0)));
	windowPoints.push_back(Point(round(x + radius * scale < WINDOW_SIZE ? x + radius * scale : WINDOW_SIZE), round(y - radius * scale > 0 ? y - radius * scale : 0)));
	windowPoints.push_back(Point(round(x - radius * scale > 0 ? x - radius * scale : 0), round(y + radius * scale < WINDOW_SIZE ? y + radius * scale : WINDOW_SIZE)));
	windowPoints.push_back(Point(round(x + radius * scale < WINDOW_SIZE ? x + radius * scale : WINDOW_SIZE), round(y + radius * scale < WINDOW_SIZE ? y + radius * scale : WINDOW_SIZE)));
	return windowPoints;
}

void Components::setX(int x) {
	this->x = x;
}

void Components::setY(int y) {
	this->y = y;
}

void Components::setType(int type) {
	this->type = type;
}

void Components::setRadius(int radius) {
	this->radius = radius;
}

void Components::setScale(double scale) {
	this->scale = scale;
}

void Components::setOrientation(double orientation) {
	this->orientation = orientation;
}

void Components::setRandFlag(bool Flag) {
	this->RandFlag = Flag;
}

void Components::setNegSampleFlag(bool Flag) {
	this->negSampleFlag = Flag;
}

int Components::getX() {
	return this->x;
}

int Components::getY() {
	return this->y;
}

int Components::getType() {
	return this->type;
}

int Components::getRadius() {
	return this->radius;
}

double Components::getScale() {
	return this->scale;
}

double Components::getOrientation() {
	return this->orientation;
}

Mat Components::drawPic() {
	return Mat();
};


// if (type == ?){...}
Mat Face::drawPic(Mat panel) {
	if (negSampleFlag == true) {
		setX(int(negSamples_x_dist(generator)));
		setY(int(negSamples_x_dist(generator)));
		setScale(face_scale_dist(generator));
		setOrientation(face_orientation_dist(generator));
	}
	else if (RandFlag == true) {
		// setX(int(face_x_dist(generator)));
		// setY(int(face_y_dist(generator)));
		// setScale(face_scale_dist(generator));
		// setOrientation(face_orientation_dist(generator));
	}
	if (type == 1) {
		// Ellipsoid
		ellipse(panel, Point(x, y), Size(int(w / 2.8 * scale), int(w / 2.4 * scale)), orientation, 0.0, 360.0, Scalar(0), thickness);
	}
	else if (type == 2) {
		// sharp face
		vector<Point> pointLists;
		pointLists.push_back(Point(x - 50, y - 50));
		pointLists.push_back(Point(x - 44, y + 34));
		pointLists.push_back(Point(x, y + 50));
		pointLists.push_back(Point(x + 44, y + 34));
		pointLists.push_back(Point(x + 50, y - 50));
		Point origin = Point(x, y);

		vector<Point> newPointLists;
		for (int i = 0; i < pointLists.size(); i++) {
			newPointLists.push_back(rescaleAndRotate(pointLists[i], Point(x, y), scale, orientation));
		}
		polylines(panel, newPointLists, true, Scalar(0), thickness, 8, 0);
	}


	return panel;
}

Mat Eye::drawPic(Mat panel) {
	if (negSampleFlag == true) {
		if (eyeType == "Right") {
			setX(int(negSamples_x_dist(generator)));
			setY(int(negSamples_x_dist(generator)));
			setScale(eye_right_scale_dist(generator));
			setOrientation(eye_right_orientation_dist(generator));
		}
		else if (eyeType == "Left") {
			setX(int(negSamples_x_dist(generator)));
			setY(int(negSamples_x_dist(generator)));
			setScale(eye_left_scale_dist(generator));
			setOrientation(eye_left_orientation_dist(generator));
		}

	} 
	else if (RandFlag == true) {
		// if (eyeType == "Right") {
		// 	setX(int(eye_right_x_dist(generator)));
		// 	setY(int(eye_right_y_dist(generator)));
		// 	setScale(eye_right_scale_dist(generator));
		// 	setOrientation(eye_right_orientation_dist(generator));
		// }
		// else if (eyeType == "Left"){
		// 	setX(int(eye_left_x_dist(generator)));
		// 	setY(int(eye_left_y_dist(generator)));
		// 	setScale(eye_left_scale_dist(generator));
		// 	setOrientation(eye_left_orientation_dist(generator));
		// }
	}

	if (type == 1) {
		// ellipsoid
		ellipse(panel, Point(x, y), Size(int(w / 12 * scale), int(w / 14.5 * scale)), orientation, 0, 360, Scalar(0),thickness);
	}
	else if (type == 2) {
		ellipse(panel, Point(x, y-10), Size(int(w / 12 * scale), int(w / 14.5 * scale)), orientation, 30, 150, Scalar(0),thickness);
	}
	else if (type == 3) {
		ellipse(panel, Point(x, y+10), Size(int(w / 12 * scale), int(w / 14.5 * scale)), orientation, 210, 310, Scalar(0),thickness);
	}
	else if (type == 4) {
		vector<Point> pointLists;
		if (eyeType == "Right") {
			pointLists.push_back(Point(x - 10, y - 10));
			pointLists.push_back(Point(x, y));
			pointLists.push_back(Point(x - 10, y + 10));
			
		}
		else if (eyeType == "Left") {
			pointLists.push_back(Point(x + 10, y - 10));
			pointLists.push_back(Point(x, y));
			pointLists.push_back(Point(x + 10, y + 10));
		}
		Point origin = Point(x, y);

		vector<Point> newPointLists;
		for (int i = 0; i < pointLists.size(); i++) {
			newPointLists.push_back(rescaleAndRotate(pointLists[i], Point(x, y), scale, orientation));
		}
		polylines(panel, newPointLists, false, Scalar(0), thickness, 8, 0);
	}
	return panel;
}

Mat Ear::drawPic(Mat panel) {
	if (negSampleFlag == true) {
		if (earType == "Right") {
			setX(int(negSamples_x_dist(generator)));
			setY(int(negSamples_x_dist(generator)));
			setScale(ear_right_scale_dist(generator));
			setOrientation(ear_right_orientation_dist(generator));
		}
		else if (earType == "Left") {
			setX(int(negSamples_x_dist(generator)));
			setY(int(negSamples_x_dist(generator)));
			setScale(ear_left_scale_dist(generator));
			setOrientation(ear_left_orientation_dist(generator));
		}

	}
	else if (RandFlag == true) {
		// if (earType == "Right") {
		// 	setX(int(ear_right_x_dist(generator)));
		// 	setY(int(ear_right_y_dist(generator)));
		// 	setScale(ear_right_scale_dist(generator));
		// 	setOrientation(ear_right_orientation_dist(generator));
		// }
		// else if (earType == "Left") {
		// 	setX(int(ear_left_x_dist(generator)));
		// 	setY(int(ear_left_y_dist(generator)));
		// 	setScale(ear_left_scale_dist(generator));
		// 	setOrientation(ear_left_orientation_dist(generator));
		// }
	}

	if (type == 1) {
		vector<Point> pointLists;
		pointLists.push_back(Point(x - 5, y - 8));
		pointLists.push_back(Point(x - 5, y + 6));
		pointLists.push_back(Point(x + 5, y + 8));
		pointLists.push_back(Point(x + 5, y - 6));

		if (earType == "Right") {
			pointLists = mirror(pointLists,	x);
		}

		vector<Point> newPointLists;
		for (int i = 0; i < pointLists.size(); i++) {
			newPointLists.push_back(rescaleAndRotate(pointLists[i], Point(x, y), scale, orientation));
		}
		polylines(panel, newPointLists, true, Scalar(0), 1, 8, 0);
	}
	else if (type == 2) {
		if (earType == "Left") {
			ellipse(panel, Point(x-7, y), Size(int(w / 12 * scale), int(w / 14.5 * scale)), orientation, -90, 90, Scalar(0), thickness);
		}
		
		if (earType == "Right") {
			ellipse(panel, Point(x+7, y), Size(int(w / 12 * scale), int(w / 14.5 * scale)), orientation, 90, 270, Scalar(0), thickness);
		}
	}
	else if (type == 3) {
		if (earType == "Left") {
			ellipse(panel, Point(x - 5, y+ int(w / 14.5 * scale)), Size(int(w / 12 * scale), int(w / 14.5 * scale)), orientation, -90, 90, Scalar(0), thickness);
			ellipse(panel, Point(x - 5, y- int(w / 14.5 * scale)), Size(int(w / 12 * scale), int(w / 14.5 * scale)), orientation, -90, 90, Scalar(0), thickness);
		}

		if (earType == "Right") {
			ellipse(panel, Point(x + 5, y + int(w / 14.5 * scale)), Size(int(w / 12 * scale), int(w / 14.5 * scale)), orientation, 90, 270, Scalar(0), thickness);
			ellipse(panel, Point(x + 5, y - int(w / 14.5 * scale)), Size(int(w / 12 * scale), int(w / 14.5 * scale)), orientation, 90, 270, Scalar(0), thickness);
		}
	}

	
	return panel;
}

Mat Nose::drawPic(Mat panel) {
	if (negSampleFlag == true) {
		setX(int(negSamples_x_dist(generator)));
		setY(int(negSamples_x_dist(generator)));
		setScale(nose_scale_dist(generator));
		setOrientation(nose_orientation_dist(generator));
	}
	else if (RandFlag == true) {
		// setX(int(nose_x_dist(generator)));
		// setY(int(nose_y_dist(generator)));
		// setScale(nose_scale_dist(generator));
		// setOrientation(nose_orientation_dist(generator));
	}

	if (type == 1) {
		// triangle
		vector<Point> pointLists;
		pointLists.push_back(Point(x, y - 16));
		pointLists.push_back(Point(x - 10, y + 10));
		pointLists.push_back(Point(x + 10, y + 10));

		vector<Point> newPointLists;
		for (int i = 0; i < pointLists.size(); i++) {
			newPointLists.push_back(rescaleAndRotate(pointLists[i], Point(x, y), scale, orientation));
		}
		polylines(panel, newPointLists, true, Scalar(0), 1, 8, 0);
	}
	else if (type == 2) {
		circle(panel, Point(x, y), 5, Scalar(0), thickness);
	}
	else if (type == 3) {
		vector<Point> pointLists;
		pointLists.push_back(Point(x -  5, y +  10));
		pointLists.push_back(Point(x -  5, y -  10));
		pointLists.push_back(Point(x +  5, y -  10));
		pointLists.push_back(Point(x +  5, y +  10));

		vector<Point> newPointLists;
		for (int i = 0; i < pointLists.size(); i++) {
			newPointLists.push_back(rescaleAndRotate(pointLists[i], Point(x, y), scale, orientation));
		}
		polylines(panel, newPointLists, true, Scalar(0), thickness, 8, 0);
	}
	return panel;
}

Mat Mouth::drawPic(Mat panel) {
	if (negSampleFlag == true) {
		setX(int(negSamples_x_dist(generator)));
		setY(int(negSamples_x_dist(generator)));
		setScale(mouth_scale_dist(generator));
		setOrientation(mouth_orientation_dist(generator));
	}
	else if (RandFlag == true) {
		// setX(int(mouth_x_dist(generator)));
		// setY(int(mouth_y_dist(generator)));
		// setScale(mouth_scale_dist(generator));
		// setOrientation(mouth_orientation_dist(generator));
	}
	if (type == 1) {
		ellipse(panel, Point(x, y), Size(int(w / 6 * scale), int(w / 6.5 * scale)), orientation, 45, 135, Scalar(0), thickness);
	}
	else if (type == 2) {
		ellipse(panel, Point(x, y+20), Size(int(w / 6 * scale), int(w / 6.5 * scale)), orientation, 210, 330, Scalar(0), thickness);
	}
	else if (type == 3) {
		vector<Point> pointLists;
		pointLists.push_back(Point(x - 10, y+5));
		pointLists.push_back(Point(x - 5, y +10));
		pointLists.push_back(Point(x, y+5));
		pointLists.push_back(Point(x + 5, y +10));
		pointLists.push_back(Point(x + 15, y+5));

		vector<Point> newPointLists;
		for (int i = 0; i < pointLists.size(); i++) {
			newPointLists.push_back(rescaleAndRotate(pointLists[i], Point(x, y), scale, orientation));
		}
		polylines(panel, newPointLists, false, Scalar(0), thickness, 8, 0);
	}
	return panel;
}
