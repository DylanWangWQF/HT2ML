//
//  DatasetOperations.h
//  HT2ML
//
//  Created by Qifan Wang on 14/01/22.
//

#ifndef DatasetOperations_h
#define DatasetOperations_h

//#include "mnist/include/mnist/mnist_reader_less.hpp"
#include "mnist/include/mnist/mnist_reader.hpp"

#include <string>
#include <vector>

using namespace std;

/*** DATASET OPERATION ***/

template <template <typename...> class Container, typename Image>
void normalize(Container<Image>& images);

vector<vector<float>> loadMnistTestImages(const string& dataset_dir, const size_t& test_limit);
vector<unsigned char> loadMnistTestLabels(const string& dataset_dir, const size_t& test_limit);

#endif /* DatasetOperations_h */
