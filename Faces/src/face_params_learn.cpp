#include <vector>
#include <random>
#include "face_params.h"

std::uniform_real_distribution<double> initial_dist_xy(20, 100);
std::default_random_engine init_generator;

double face_x_mean = initial_dist_xy(init_generator);
double face_y_mean = initial_dist_xy(init_generator);
double nose_x_mean = initial_dist_xy(init_generator);
double nose_y_mean = initial_dist_xy(init_generator);
double mouth_x_mean = initial_dist_xy(init_generator);
double mouth_y_mean = initial_dist_xy(init_generator);
double ear_left_x_mean = initial_dist_xy(init_generator);
double ear_left_y_mean = initial_dist_xy(init_generator);
double ear_right_x_mean = initial_dist_xy(init_generator);
double ear_right_y_mean = initial_dist_xy(init_generator);
double eye_left_x_mean = initial_dist_xy(init_generator);
double eye_left_y_mean = initial_dist_xy(init_generator);
double eye_right_x_mean = initial_dist_xy(init_generator);
double eye_right_y_mean = initial_dist_xy(init_generator);

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