//
// Created by kore on 9/8/21.
//

#include "readData.h"

int reverse_int(int i){
    unsigned char ch1, ch2, ch3, ch4;
    ch1 = i & 255;
    ch2 = (i >> 8) & 255;
    ch3 = (i >> 16) & 255;
    ch4 = (i >> 24) & 255;
    return ((int) ch1 << 24) + ((int) ch2 << 16) + ((int) ch3 << 8) + ch4;
}



//template <typename T>
void read_MNIST_data(bool train, std::vector<std::vector<double> > &vec, int& number_of_images, int& number_of_features, int NUM_IMAGES){
    std::ifstream file;

    double check1, check2;

    if (train == true)
        file.open("../Datasets/train-images-idx3-ubyte",std::ios::binary);
    else
        file.open("../Datasets/t10k-images-idx3-ubyte",std::ios::binary);


    if(!file){
        std::cout<<"Unable to open file";
        return;
    }
    if (file.is_open()){
        int magic_number = 0;
        int n_rows = 0;
        int n_cols = 0;
        file.read((char*) &magic_number, sizeof(magic_number));
        magic_number = reverse_int(magic_number);
        file.read((char*) &number_of_images, sizeof(number_of_images));
        number_of_images = reverse_int(number_of_images)/10;

        //number_of_images = (int) (number_of_images /5);

        if(train == true)
            number_of_images = NUM_IMAGES;
        file.read((char*) &n_rows, sizeof(n_rows));
        n_rows = reverse_int(n_rows);
        file.read((char*) &n_cols, sizeof(n_cols));
        n_cols = reverse_int(n_cols);
        number_of_features = n_rows * n_cols;
        std::cout << "Number of Images: " << number_of_images << std::endl;
        std::cout << "Number of Features: " << number_of_features << std::endl;
        int ft=1;
        for(int i = 0; i < number_of_images; ++i){
            std::vector<double> tp;
            for(int r = 0; r < n_rows; ++r)
                for(int c = 0; c < n_cols; ++c){
                    unsigned char temp = 0;
                    file.read((char*) &temp, sizeof(temp));
                    tp.push_back(((double) temp));
/*
                    check1 = (double) temp;
                    check2 = (double) ((int) temp);

                    if(check1 != check2)cout<< check1 << endl;*/
                    //if(ft){
                    //    cout<< check1 << endl;

                    //}
                }
            //ft = 0;
            vec.push_back(tp);
        }
    }
}

//template <typename T>
void read_MNIST_labels(bool train, std::vector<double> &vec, int NUM_IMAGES){
    std::ifstream file;
    if (train == true)
        file.open("../Datasets/train-labels-idx1-ubyte",std::ios::binary);
    else
        file.open("../Datasets/t10k-labels-idx1-ubyte", std::ios::binary);
    if(!file){
        std::cout << "Unable to open file";
        return;
    }
    if (file.is_open()){
        int magic_number = 0;
        int number_of_images = 0;
        file.read((char*) &magic_number, sizeof(magic_number));
        magic_number = reverse_int(magic_number);
        file.read((char*) &number_of_images, sizeof(number_of_images));
        number_of_images = reverse_int(number_of_images)/10;
        if(train == true)
            number_of_images = NUM_IMAGES;
        for(int i = 0; i < number_of_images; ++i){
            unsigned char temp = 0;
            file.read((char*) &temp, sizeof(temp));
            //vec.push_back((T) temp);
            //cout << (int) temp <<endl;

            if((double) temp == 2)
                vec.push_back((double) 0);
                //else if((T) temp == 9)
                //    vec.push_back((T) 1);
            else
                vec.push_back((double) 1);
        }
    }
}

