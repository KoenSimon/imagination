//
// Created by kore on 9/8/21.
//

#ifndef PPKNN_READDATA_H
#define PPKNN_READDATA_H

#include <iostream>
#include <vector>
#include <fstream>


int reverse_int(int i);
void read_MNIST_data(bool train, std::vector<std::vector<double> > &vec, int& number_of_images, int& number_of_features, int NUM_IMAGES);
void read_MNIST_labels(bool train, std::vector<double> &vec, int NUM_IMAGES);

#endif //PPKNN_READDATA_H
