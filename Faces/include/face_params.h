#ifndef FACE_PARAMS_H
#define FACE_PARAMS_H

#include <vector>

// Position mean
extern std::vector<double> param_reference;
extern double face_x_mean;
extern double face_y_mean;
extern double nose_x_mean;
extern double nose_y_mean;
extern double mouth_x_mean;
extern double mouth_y_mean;
extern double ear_left_x_mean;
extern double ear_left_y_mean;
extern double ear_right_x_mean;
extern double ear_right_y_mean;
extern double eye_left_x_mean;
extern double eye_left_y_mean;
extern double eye_right_x_mean;
extern double eye_right_y_mean;

// Position variance
extern double face_x_var;
extern double face_y_var;
extern double nose_x_var;
extern double nose_y_var;
extern double mouth_x_var;
extern double mouth_y_var;
extern double ear_left_x_var;
extern double ear_left_y_var;
extern double ear_right_x_var;
extern double ear_right_y_var;
extern double eye_left_x_var;
extern double eye_left_y_var;
extern double eye_right_x_var;
extern double eye_right_y_var;

extern int face_x_mean_idx;
extern int face_y_mean_idx;
extern int nose_x_mean_idx;
extern int nose_y_mean_idx;
extern int mouth_x_mean_idx;
extern int mouth_y_mean_idx;
extern int ear_left_x_mean_idx;
extern int ear_left_y_mean_idx;
extern int ear_right_x_mean_idx;
extern int ear_right_y_mean_idx;
extern int eye_left_x_mean_idx;
extern int eye_left_y_mean_idx;
extern int eye_right_x_mean_idx;
extern int eye_right_y_mean_idx;

extern int face_x_var_idx;
extern int face_y_var_idx;
extern int nose_x_var_idx;
extern int nose_y_var_idx;
extern int mouth_x_var_idx;
extern int mouth_y_var_idx;
extern int ear_left_x_var_idx;
extern int ear_left_y_var_idx;
extern int ear_right_x_var_idx;
extern int ear_right_y_var_idx;
extern int eye_left_x_var_idx;
extern int eye_left_y_var_idx;
extern int eye_right_x_var_idx;
extern int eye_right_y_var_idx;

// Scale variance
extern double face_scale_var;
extern double nose_scale_var;
extern double mouth_scale_var;
extern double ear_left_scale_var;
extern double ear_right_scale_var;
extern double eye_left_scale_var;
extern double eye_right_scale_var;

// Orientation variance
extern double face_pose_var;
extern double nose_pose_var;
extern double mouth_pose_var;
extern double ear_left_pose_var;
extern double ear_right_pose_var;
extern double eye_left_pose_var;
extern double eye_right_pose_var;


// radius
extern int face_radius;
extern int eye_left_radius;
extern int eye_right_radius;
extern int ear_left_radius;
extern int ear_right_radius;
extern int nose_radius;
extern int mouth_radius;

// Negtive Samples range
extern int lower_limit;
extern int higher_limit;

#endif // FACE_PARAMS_H
