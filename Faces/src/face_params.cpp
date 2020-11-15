#include <vector>
#include "face_params.h"

double face_x_mean = 64;
double face_y_mean = 64;
double nose_x_mean = 64;
double nose_y_mean = 70;
double mouth_x_mean = 64;
double mouth_y_mean = 85;
double ear_left_x_mean = 118;
double ear_left_y_mean = 50;
double ear_right_x_mean = 10;
double ear_right_y_mean = 50;
double eye_left_x_mean = 80;
double eye_left_y_mean = 42;
double eye_right_x_mean = 48;
double eye_right_y_mean = 42;

std::vector<double> param_reference =
 {face_x_mean, face_y_mean, nose_x_mean, nose_y_mean, mouth_x_mean, mouth_y_mean,
  ear_left_x_mean, ear_left_y_mean, ear_right_x_mean, ear_right_y_mean, eye_left_x_mean,
  eye_left_y_mean, eye_right_x_mean, eye_right_y_mean};

int face_x_mean_idx = 0;
int face_y_mean_idx = 1;
int nose_x_mean_idx = 2;
int nose_y_mean_idx = 3;
int mouth_x_mean_idx = 4;
int mouth_y_mean_idx = 5;
int ear_left_x_mean_idx = 6;
int ear_left_y_mean_idx = 7;
int ear_right_x_mean_idx = 8;
int ear_right_y_mean_idx = 9;
int eye_left_x_mean_idx = 10;
int eye_left_y_mean_idx = 11;
int eye_right_x_mean_idx = 12;
int eye_right_y_mean_idx = 13;